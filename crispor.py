#!/usr/bin/env python2.7
# the tefor crispr tool
# can be run as a CGI or from the command line

# OOF scores are WRONG for Cpf1! -> where is the cut site?

# python std library
import subprocess, tempfile, logging, atexit, glob, shutil
import Cookie, time, sys, cgi, re, random, platform, os
import hashlib, base64, string, operator, urllib, sqlite3, time
import traceback, json, pwd, pickle
import common_functions
import parser_functions
import get_offtargets

from bio_functions import *
from fasta_functions import *
from common_functions import *
from basic_imports import *
from constants import *
from collections import defaultdict, namedtuple
from datetime import datetime
from itertools import product
from os.path import join, isfile,  dirname, isdir, abspath
from StringIO import StringIO
 




# ====== FUNCTIONS =====
contentLineDone = False

def rulerString(maxLen):
    " return line with positions every 10 chars "
    texts = []
    for i in range(0, maxLen, 10):
        numStr = str(i)
        texts.append(numStr)
        spacer = "".join([" "]*(10-len(numStr)))
        texts.append(spacer)
    return "".join(texts)



def iterOneDelSeqs(seq):
    """ given a seq, create versions with each bp removed. Avoid duplicates 
    yields (delPos, seq)
    >>> list(iterOneDelSeqs("AATGG"))
    [(0, 'ATGG'), (2, 'AAGG'), (3, 'AATG')]
    """
    doneSeqs = set()
    for i in range(0, len(seq)):
        delSeq = seq[:i]+seq[i+1:]
        if delSeq not in doneSeqs:
            yield i, delSeq
        doneSeqs.add(delSeq)

def flankSeqIter(seq, startDict, pam, doFilterNs):
    """ given a seq and dictionary of pos -> strand and the length of the pamSite
    yield tuples of (name, pamStart, guideStart, strand, flankSeq, pamSeq)

    if doFilterNs is set, will not return any sequences that contain an N character
    """
    pamLen = len(pam)
    GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(pam)
    startList = sorted(startDict.keys())
    for pamStart in startList:
        strand = startDict[pamStart]

        if cpf1Mode: # Cpf1: get the sequence to the right of the PAM
            if strand=="+":
                guideStart = pamStart+pamLen
                flankSeq = seq[guideStart:guideStart+GUIDELEN]
                pamSeq = seq[pamStart:pamStart+pamLen]
            else: # strand is minus
                guideStart = pamStart-GUIDELEN
                flankSeq = revComp(seq[guideStart:pamStart])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])
        else: # common case: get the sequence on the left side of the PAM
            if strand=="+":
                guideStart = pamStart-GUIDELEN
                flankSeq = seq[guideStart:pamStart]
                pamSeq = seq[pamStart:pamStart+pamLen]
            else: # strand is minus
                guideStart = pamStart+pamLen
                flankSeq = revComp(seq[guideStart:guideStart+GUIDELEN])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])

        if "N" in flankSeq and doFilterNs:
            continue

        yield "s%d%s" % (pamStart, strand), pamStart, guideStart, strand, flankSeq, pamSeq


def highlightMismatches(guide, offTarget, pamLen):
    " return a string that marks mismatches between guide and offtarget with * "
    if cpf1Mode:
        offTarget = offTarget[pamLen:]
    else:
        offTarget = offTarget[:-pamLen]
    assert(len(guide)==len(offTarget))

    s = []
    for x, y in zip(guide, offTarget):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)

def makeAlnStr(seq1, seq2, pam, mitScore, cfdScore, posStr):
    " given two strings of equal length, return a html-formatted string that highlights the differences "
    lines = [ [], [], [] ]
    last12MmCount = 0

    if cpf1Mode:
        lines[0].append("<i>"+seq1[:len(pam)]+"</i> ")
        lines[1].append("<i>"+seq2[:len(pam)]+"</i> ")
        lines[2].append("".join([" "]*(len(pam)+1)))

    if cpf1Mode:
        guideStart = len(pam)
        guideEnd = len(seq1)
    else:
        guideStart = 0
        guideEnd = len(seq1)-len(pam)

    for i in range(guideStart, guideEnd):
        if seq1[i]==seq2[i]:
            lines[0].append(seq1[i])
            lines[1].append(seq2[i])
            lines[2].append(" ")
        else:
            lines[0].append("<b>%s</b>"%seq1[i])
            lines[1].append("<b>%s</b>"%seq2[i])
            lines[2].append("*")
            if i>7:
                last12MmCount += 1

    if not cpf1Mode:
        lines[0].append(" <i>"+seq1[-len(pam):]+"</i>")
        lines[1].append(" <i>"+seq2[-len(pam):]+"</i>")
    lines = ["".join(l) for l in lines]

    if len(posStr)>1 and posStr[0].isdigit():
        posStr = "chr"+posStr

    htmlText1 = "<small><pre>guide:      %s<br>off-target: %s<br>            %s</pre>" \
        % (lines[0], lines[1], lines[2])
    if cpf1Mode:
        htmlText2 = "CPf1: No off-target scores available</small>"
    else:
        if cfdScore==None:
            cfdStr = "Cannot calculate CFD score on non-ACTG characters"
        else:
            cfdStr = "%f" % cfdScore
        htmlText2 = "CFD Off-target score: %s<br>MIT Off-target score: %.2f<br>Position: %s</small>" % (cfdStr, mitScore, posStr)
    hasLast12Mm = last12MmCount>0
    return htmlText1+htmlText2, hasLast12Mm

def parsePos(text):
    """ parse a string of format chr:start-end:strand and return a 4-tuple
    Strand defaults to + and end defaults to start+23
    """
    if text!=None and len(text)!=0 and text!="?":
        fields = text.split(":")
        if len(fields)==2:
            chrom, posRange = fields
            strand = "+"
        else:
            chrom, posRange, strand = fields
        posRange = posRange.replace(",","")
        start, end = posRange.split("-")
        start, end = int(start), int(end)
    else:
        chrom, start, end, strand = "", 0, 0, "+"
    return chrom, start, end, strand

def makePosList(countDict, guideSeq, pam, inputPos):
    """ for a given guide sequence, return a list of tuples that
    describes the offtargets sorted by score and a string to describe the offtargets in the
    format x/y/z/w of mismatch counts
    inputPos has format "chrom:start-end:strand". All 0MM matches in this range
    are ignored from scoring ("ontargets")
    Also return the same description for just the last 12 bp and the score
    of the guide sequence (calculated using all offtargets).
    """
    inChrom, inStart, inEnd, inStrand = parsePos(inputPos)
    count = 0
    otCounts = []
    posList = []
    mitOtScores = []
    cfdScores = []
    last12MmCounts = []
    ontargetDesc = ""
    subOptMatchCount = 0

    # for each edit distance, get the off targets and iterate over them
    foundOneOntarget = False
    for editDist in range(0, maxMMs+1):
        #print countDict,"<p>"
        matches = countDict.get(editDist, [])

        #print otCounts,"<p>"
        last12MmOtCount = 0

        # create html and score for every offtarget
        otCount = 0
        for chrom, start, end, otSeq, strand, segType, geneNameStr, x1Count in matches:
            # skip on-targets
            if segType!="":
                segTypeDesc = segTypeConv[segType]
                geneDesc = segTypeDesc+":"+geneNameStr
                geneDesc = geneDesc.replace("|", "-")
            else:
                geneDesc = geneNameStr

            # is this the on-target?
            # if we got a genome position, use it. Otherwise use a random off-target with 0MMs
            # as the on-target
            if editDist==0 and \
                ((chrom==inChrom and start >= inStart and end <= inEnd and x1Count < MAXOCC) \
                or (inChrom=='' and foundOneOntarget==False and x1Count < MAXOCC)):
                foundOneOntarget = True
                ontargetDesc = geneDesc
                continue

            otCount += 1
            guideNoPam = guideSeq[:len(guideSeq)-len(pam)]
            otSeqNoPam = otSeq[:len(otSeq)-len(pam)]
            if len(otSeqNoPam)==19:
                otSeqNoPam = "A"+otSeqNoPam # should not change the score a lot, weight0 is very low
            if pam!="TTTN":
                # MIT score must not include the PAM
                mitScore = calcHitScore(guideNoPam, otSeqNoPam)
                # CFD score must include the PAM
                cfdScore = calcCfdScore(guideSeq, otSeq)
            else:
                mitScore=0.0
                cfdScore=0.0

            mitOtScores.append(mitScore)
            if cfdScore != None:
                cfdScores.append(cfdScore)

            posStr = "%s:%d-%s:%s" % (chrom, int(start)+1,end, strand)
            alnHtml, hasLast12Mm = makeAlnStr(guideSeq, otSeq, pam, mitScore, cfdScore, posStr)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (otSeq, mitScore, cfdScore, editDist, posStr, geneDesc, alnHtml) )
            # taking the maximum is probably not necessary, 
            # there should be only one offtarget for X1-exceeding matches
            subOptMatchCount = max(int(x1Count), subOptMatchCount)

        last12MmCounts.append(str(last12MmOtCount))
        # create a list of number of offtargets for this edit dist
        otCounts.append( str(otCount) )

    # calculate the guide scores
    if pam=="TTTN":
        guideScore = -1
        guideCfdScore = -1
    else:
        if subOptMatchCount > MAXOCC:
            guideScore = 0
            guideCfdScore = 0
        else:
            guideScore = calcMitGuideScore(sum(mitOtScores))
            guideCfdScore = calcMitGuideScore(sum(cfdScores))

    # obtain the off-target info: coordinates, descriptions and off-target counts
    if subOptMatchCount > MAXOCC:
        posList = []
        ontargetDesc = ""
        last12DescStr = ""
        otDescStr = ""
    else:
        otDescStr = "&thinsp;-&thinsp;".join(otCounts)
        last12DescStr = "&thinsp;-&thinsp;".join(last12MmCounts)

    if pam=="TTTN":
        # sort by edit dist if using Cfp1
        posList.sort(key=operator.itemgetter(3))
    else:
        # sort by CFD score if we have it
        posList.sort(reverse=True, key=operator.itemgetter(2))

    return posList, otDescStr, guideScore, guideCfdScore, last12DescStr, \
        ontargetDesc, subOptMatchCount

# --- START OF SCORING ROUTINES 

# MIT offtarget scoring

# aka Matrix "M"
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    if len(string1)==21 and len(string2)==21:
        string1 = string1[-20:]
        string2 = string2[-20:]

    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about 
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====
# The CFD score is an improved specificity score 

def get_mm_pam_scores():
    """
    """
    dataDir = join(dirname(__file__), 'CFD_Scoring')
    mm_scores = pickle.load(open(join(dataDir, 'mismatch_score.pkl'),'rb'))
    pam_scores = pickle.load(open(join(dataDir, 'pam_scores.pkl'),'rb'))
    return (mm_scores,pam_scores)

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

mm_scores, pam_scores = None, None

def calcCfdScore(guideSeq, otSeq):
    """ based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    """
    global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores,pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:20]
        cfd_score = calc_cfd(wt,sg,pam)
        return cfd_score
# ==== END CFD score source provided by John Doench

# --- END OF SCORING ROUTINES 

def parseChromSizes(genome,genomesDir):
    " return chrom sizes as dict chrom -> size "
    genomeDir = genomesDir # make local
    sizeFname = "%(genomeDir)s/%(genome)s/%(genome)s.sizes" % locals()
    ret = {}
    for line in open(sizeFname).read().splitlines():
        fields = line.split()
        chrom, size = fields[:2]
        ret[chrom] = int(size)
    return ret

def extendAndGetSeq(db,genomesDir, chrom, start, end, strand, flank=100):
    """ extend (start, end) by flank and get sequence for it using twoBitTwoFa.
    Return None if not possible to extend.
    #>>> extendAndGetSeq("hg19", "chr21", 10000000, 10000005, "+", flank=3)
    #'AAGGAATGTAG'
    """
    chromSizes = parseChromSizes(db,genomesDir)
    maxEnd = chromSizes[chrom]+1

    start -= flank
    end += flank
    if start < 0 or end > maxEnd:
        return None

    genomeDir = genomesDir
    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    progDir = binDir
    genome = db
    cmd = "%(progDir)s/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -seq=%(chrom)s -start=%(start)s -end=%(end)s" % locals()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()
    faFile = StringIO(seqStr)
    seqs = parseFasta(faFile)
    assert(len(seqs)==1)
    seq = seqs.values()[0].upper()
    if strand=="-":
        seq = revComp(seq)
    return seq

def getExtSeq(seq, start, end, strand, extUpstream, extDownstream, extSeq=None, extFlank=100):
    """ extend (start,end) by extUpstream and extDownstream and return the subsequence
    at this position in seq.
    Return None if there is not enough space to extend (start, end).
    extSeq is a sequence with extFlank additional flanking bases on each side. It can be provided
    optionally and is used if needed to return a subseq.
    Careful: returned sequence might contain lowercase letters.
    >>> getExtSeq("AACCTTGG", 2, 4, "+", 2, 4)
    'AACCTTGG'
    >>> getExtSeq("CCAACCTTGGCC", 4, 6, "-", 2, 3)
    'AAGGTTG'
    >>> getExtSeq("AA", 0, 2, "+", 2, 3)
    >>> getExtSeq("AA", 0, 2, "+", 2, 3, extSeq="CAGAATGA", extFlank=3)
    'AGAATGA'
    >>> getExtSeq("AA", 0, 2, "-", 2, 3, extSeq="CAGAATGA", extFlank=3)
    'CATTCTG'
    """
    assert(start>=0)
    assert(end<=len(seq))
    # check if the extended sequence really contains the whole input seq 
    # e.g. when user has added nucleotides to a otherwise matching sequence
    if extSeq!=None and (seq.upper() not in extSeq.upper()):
        debug("seq is not in extSeq")
        extSeq = None

    # extend
    if strand=="+":
        extStart, extEnd = start-extUpstream, end+extDownstream
    else:
        extStart, extEnd = start-extDownstream, end+extUpstream

    # check for out of bounds and get seq
    if extStart >= 0 and extEnd <= len(seq):
        subSeq = seq[extStart:extEnd]
    else:
        if extSeq==None:
            return None
        # lift to extSeq coords and get seq
        extStart += extFlank
        extEnd += extFlank
        assert(extStart >= 0)
        assert(extEnd <= len(extSeq))
        subSeq = extSeq[extStart:extEnd]

    if strand=="-":
        subSeq = revComp(subSeq)

    return subSeq

def pamStartToGuideRange(startPos, strand, pamLen,cpf1Mode,GUIDELEN):
    """ given a PAM start position and its strand, return the (start,end) of the guide.
    Coords can be negative or exceed the length of the input sequence.
    """
    if not cpf1Mode:
        if strand=="+":
            return (startPos-GUIDELEN, startPos)
        else: # strand is minus
            return (startPos+pamLen, startPos+pamLen+GUIDELEN)
    else:
        if strand=="+":
            return (startPos+pamLen, startPos+pamLen+GUIDELEN)
        else: # strand is minus
            return (startPos-GUIDELEN, startPos)



def readEnzymes():
    """ parse restrSites.txt and
    return as dict length -> list of (name, suppliers, seq) """
    fname = "restrSites2.txt"
    enzList = {}
    for line in open(join(baseDir, fname)):
        if line.startswith("#"):
            continue
        seq, name, suppliers = line.rstrip("\n").split("\t")
        suppliers = tuple(suppliers.split(","))
        enzList.setdefault(len(seq), []).append( (name, suppliers, seq) )
    return enzList
        
def patMatch(seq, pat, notDegPos=None):
    """ return true if pat matches seq, both have to be same length 
    do not match degenerate codes at position notDegPos (0-based)
    """
    assert(len(seq)==len(pat))
    for x in range(0, len(pat)):
        patChar = pat[x]
        nuc = seq[x]

        assert(patChar in "MKYRACTGNWSD")
        assert(nuc in "MKYRACTGNWSD")

        if notDegPos!=None and x==notDegPos and patChar!=nuc:
            #print x, seq, pat, notDegPos, patChar, nuc, "<br>"
            return False

        if patChar=="N":
            continue
        if patChar=="D" and nuc in ["AGT"]:
            continue
        if patChar=="W" and nuc in ["A", "T"]:
            continue
        if patChar=="S" and nuc in ["G", "C"]:
            continue
        if patChar=="M" and nuc in ["A", "C"]:
            continue
        if patChar=="K" and nuc in ["T", "G"]:
            continue
        if patChar=="R" and nuc in ["A", "G"]:
            continue
        if patChar=="Y" and nuc in ["C", "T"]:
            continue
        if patChar!=nuc:
            return False
    return True

def findSite(seq, restrSite):
    """ return the positions where restrSite matches seq 
    seq can be longer than restrSite
    Do not allow degenerate characters to match at position len(restrSite) in seq
    """
    posList = []
    for i in range(0, len(seq)-len(restrSite)+1):
        subseq = seq[i:i+len(restrSite)]
        #print subseq==restrSite, subseq, restrSite,"<br>"

        # JP does not want any potential site to be suppressed
        #if i<len(restrSite):
            #isMatch = patMatch(subseq, restrSite, len(restrSite)-i-1)
        #else:
            #isMatch = patMatch(subseq, restrSite)
        isMatch = patMatch(subseq, restrSite)

        if isMatch:
            posList.append( (i, i+len(restrSite)) )
    return posList

def matchRestrEnz(allEnzymes, guideSeq, pamSeq):
    """ return list of enzymes that overlap the -3 position in guideSeq
    returns dict (name, pattern, suppliers) -> list of matching positions
    """
    matches = defaultdict(set)
    #print guideSeq, pamSeq, "<br>"
    fullSeq = concatGuideAndPam(guideSeq, pamSeq)

    for siteLen, sites in allEnzymes.iteritems():
        if cpf1Mode:
            # most modified position: 4nt from the end
            # see http://www.nature.com/nbt/journal/v34/n8/full/nbt.3620.html
            # Figure 1
            startSeq = len(fullSeq)-4-(siteLen)+1
        else:
            # most modified position for Cas9: 3bp from the end
            startSeq = len(fullSeq)-len(pamSeq)-3-(siteLen)+1

        seq = fullSeq[startSeq:]
        for name, suppliers, restrSite in sites:
            posList = findSite(seq, restrSite)
            if len(posList)!=0:
                liftOffset = startSeq
                posList = [(liftOffset+x, liftOffset+y) for x,y in posList]
                matches.setdefault((name, restrSite, suppliers), set()).update(posList)
    return matches

def mergeGuideInfo(seq, startDict, pamPat, otMatches, inputPos, effScores, sortBy=None):
    """
    merges guide information from the sequence, the efficiency scores and the off-targets.
    creates rows with fields:


    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "effScore", "mhScore", "oofScore" or "pos"
    """
    allEnzymes = readEnzymes()

    guideData = []
    guideScores = {}
    hasNotFound = False

    pamSeqs = list(flankSeqIter(seq, startDict, pamPat, True))

    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = concatGuideAndPam(guideSeq, pamSeq)
            mutEnzymes = matchRestrEnz(allEnzymes, guideSeq, pamSeq)
            posList, otDesc, guideScore, guideCfdScore, last12Desc, ontargetDesc, \
               subOptMatchCount = \
                   makePosList(pamMatches, guideSeqFull, pamPat, inputPos)

        # no off-targets found?
        else:
            posList, otDesc, guideScore = None, "Not found", None
            guideCfdScore = None
            last12Desc = ""
            hasNotFound = True
            mutEnzymes = []
            ontargetDesc = ""
            subOptMatchCount = False
            seq34Mer = None

        guideRow = [guideScore, guideCfdScore, effScores.get(pamId, {}), pamStart, guideStart, strand, pamId, guideSeq, pamSeq, posList, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount]
        guideData.append( guideRow )
        guideScores[pamId] = guideScore

    if sortBy == "pos":
        sortFunc = (lambda row: row[3])
        reverse = False
    elif sortBy is not None and sortBy!="spec":
        sortFunc = (lambda row: row[2].get(sortBy, 0))
        reverse = True
    else:
        sortFunc = operator.itemgetter(0)
        reverse = True

    guideData.sort(reverse=reverse, key=sortFunc)

    return guideData, guideScores, hasNotFound

def hasGeneModels(org):
    " return true if this organism has gene model information "
    geneFname = join(genomesDir, org, org+".segments.bed")
    return isfile(geneFname)

def scoreToColor(guideScore):
    if guideScore > 50:
        color = "#32cd32"
    elif guideScore > 20:
        color = "#ffff00"
    elif guideScore==-1:
        color = "black"
    else:
        color = "#aa0114"
    return color


def filterOts(otDatas, minScore):
    " remove all offtargets with score < minScore "
    newList = []
    for otData in otDatas:
        score = otData[1]
        if score > minScore:
            newList.append(otData)
    return newList

def findOtCutoff(otData):
    " try cutoffs 0.5, 1.0, 2.0, 3.0 until not more than 20 offtargets left "
    for cutoff in [0.3, 0.5, 1.0, 2.0, 3.0, 10.0, 99.9]:
        otData = filterOts(otData, cutoff)
        if len(otData)<=30:
            return otData, cutoff

    if len(otData)>30:
        return otData[:30], None

    return otData, 1000



def firstFreeLine(lineMasks, y, start, end):
    " recursively search for first free line to place a feature (start, end) "
    #print "called with y", y
    if y>=len(lineMasks):
        return None
    lineMask = lineMasks[y]
    for x in range(start, end):
        if lineMask[x]!=0:
            return firstFreeLine(lineMasks, y+1, start, end)
        else:
            return y
    return None

def distrOnLines(seq, startDict, featLen):
    """ given a dict with start -> (start,end,name,strand) and a motif len, create lines of annotations such that
        the motifs don't overlap on the lines 
    """
    # max number of lines in y direction to draw
    MAXLINES = 18
    # amount of free space around each feature
    SLOP = 2

    # bitmask, one per line, 1 = we have a feature here, 0 = no feature here
    lineMasks = []
    for i in range(0, MAXLINES):
        lineMasks.append( [0]* (len(seq)+10) )

    # dict with lineCount (0...MAXLINES) -> list of (start, strand) tuples
    ftsByLine = defaultdict(list)
    maxY = 0
    for start in sorted(startDict):
        end = start+featLen
        strand = startDict[start]

        # Cannot use Unicode here: these symbols are not part of the
        # monospace font on some platforms and therefore their width
        # is not the same as the other characters
        #arrNE = u'\u2197'
        #arrSE = u'\u2198'
        arrNE = u'/'
        arrSE = u'\\'
        #arrNE = u'\u2a3c' # hebrew
        #arrSE = u'\ufb27' # math
        ftSeq = seq[start:end]
        if strand=="+":
            if cpf1Mode:
                label = '%s'%(ftSeq)+u'.................%s....%s' % (arrNE, arrSE)
                startFt = start
                endFt = start+len(label)
            else:
                #label = '%s..%s'%(seq[start-3].lower(), ftSeq)
                label = '---%s'%(ftSeq)
                startFt = start - 3
                endFt = end
        else:
            if cpf1Mode:
                spc1 = "...."
                spc2 = "................."
                labelPrefix = u'%s%s%s%s' % (arrSE, spc1, arrNE, spc2)
                label = labelPrefix + ftSeq
                startFt = start - len(labelPrefix)
                endFt = startFt+len(label)
            else:
                #label = '%s..%s'%(ftSeq, seq[end+2].lower())
                label = '%s---'%(ftSeq)
                startFt = start
                endFt = end + 3

        y = firstFreeLine(lineMasks, 0, startFt, endFt)
        if y==None:
            errAbort("not enough space to plot features")

        # fill the current mask
        mask = lineMasks[y]
        for i in range(max(startFt-SLOP, 0), min(endFt+SLOP, len(seq))):
            mask[i]=1

        maxY = max(y, maxY)

        pamId = "s%d%s" % (start, strand)
        ft = (startFt, endFt, label, strand, pamId) 
        ftsByLine[y].append(ft )
    return ftsByLine, maxY

def writePamFlank(seq, startDict, pam, faFname):
    " write pam flanking sequences to fasta file, optionally with versions where each nucl is removed "
    #print "writing pams to %s<br>" % faFname
    faFh = open(faFname, "w")
    for pamId, pamStart, guideStart, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, pam, True):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()


def debug(msg):
    logging.debug(msg)
    if DEBUG:
        print msg

def parseOfftargets(bedFname):
    """ parse a bed file with annotataed off target matches from overlapSelect,
    has two name fields, one with the pam position/strand and one with the
    overlapped segment 
    
    return as dict pamId -> editDist -> (chrom, start, end, seq, strand, segType, segName, x1Score)
    segType is "ex" "int" or "ig" (=intergenic)
    if intergenic, geneNameStr is two genes, split by |
    """
    # example input:
    # chrIV 9864393 9864410 s41-|-|5|ACTTGACTG|0    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5|ACTGTAGCTAGCT|9999    chrIV   9864408 9864470 in:K07F5.16
    debug("reading offtargets from %s" % bedFname)

    # first sort into dict (pamId,chrom,start,end,editDist,strand) 
    # -> (segType, segName) 
    pamData = {}
    for line in open(bedFname):
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name, segment = fields
        # hg38: ignore alternate chromosomes otherwise the 
        # regions on the main chroms look as if they could not be 
        # targeted at all with Cas9
        if chrom.endswith("_alt"):
            continue
        nameFields = name.split("|")
        pamId, strand, editDist, seq = nameFields[:4]

        if len(nameFields)>4:
            x1Count = int(nameFields[4])
        else:
            x1Count = 0
        editDist = int(editDist)
        # some gene models include colons
        if ":" in segment:
            segType, segName = string.split(segment, ":", maxsplit=1)
        else:
            segType, segName = "", segment
        start, end = int(start), int(end)
        otKey = (pamId, chrom, start, end, editDist, seq, strand, x1Count)

        # if a offtarget overlaps an intron/exon or ig/exon boundary it will
        # appear twice; in this case, we only keep the exon offtarget
        if otKey in pamData and segType!="ex":
            continue
        pamData[otKey] = (segType, segName)

    # index by pamId and edit distance
    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.iteritems():
        pamId, chrom, start, end, editDist, seq, strand, x1Score = otKey
        segType, segName = otVal
        otTuple = (chrom, start, end, seq, strand, segType, segName, x1Score)
        indexedOts[pamId].setdefault(editDist, []).append( otTuple )

    return indexedOts

class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))

def annotateBedWithPos(inBed, outBed):
    """
    given an input bed4 and an output bed filename, add an additional column 5 to the bed file
    that is a descriptive text of the chromosome pos (e.g. chr1:1.23 Mbp).
    """
    ofh = open(outBed, "w")
    for line in open(inBed):
        chrom, start = line.split("\t")[:2]
        start = int(start)

        if start>1000000:
            startStr = "%.2f Mbp" % (float(start)/1000000)
        else:
            startStr = "%.2f Kbp" % (float(start)/1000)
        desc = "%s %s" % (chrom, startStr)

        ofh.write(line.rstrip("\n"))
        ofh.write("\t")
        ofh.write(desc)
        ofh.write("\n")
    ofh.close()

def calcGuideEffScores(seq, extSeq, pam):
    """ given a sequence and an extended sequence, get all potential guides
    with pam, extend them to 100mers and score them with various eff. scores. 
    Return a
    list of rows [headers, (guideSeq, 100mer, score1, score2, score3,...), ... ]

    extSeq can be None, if we were unable to extend the sequence
    """
    GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(pam)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    startDict, endSet = findAllPams(seq, pam)
    pamInfo = list(flankSeqIter(seq, startDict, pam, False))

    guideIds = []
    guides = []
    longSeqs = []
    for pamId, startPos, guideStart, strand, guideSeq, pamSeq in pamInfo:
        guides.append(guideSeq+pamSeq)
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam),cpf1Mode,GUIDELEN)
        longSeq = getExtSeq(seq, gStart, gEnd, strand, 50-GUIDELEN, 50, extSeq)
        if longSeq!=None:
            longSeqs.append(longSeq)
            guideIds.append(pamId)

    if len(longSeqs)>0:
        if cpf1Mode:
            mh, oof, mhSeqs = crisporEffScores.calcAllBaeScores(crisporEffScores.trimSeqs(longSeqs, -50, 50))
            effScores = {}
            effScores["oof"] = oof
            effScores["mh"] = mh

        else:
            effScores = crisporEffScores.calcAllScores(longSeqs)

    else:
        effScores = {}
    scoreNames = effScores.keys()

    # reformat to rows, write all scores to file
    rows = []
    for i, (guideId, guide, longSeq) in enumerate(zip(guideIds, guides, longSeqs)):
        row = [guideId, guide, longSeq]
        for scoreName in scoreNames:
            scoreList = effScores[scoreName]
            if len(scoreList) > 0:
                row.append(scoreList[i])
            else:
                row.append("noScore?")
        rows.append(row)

    headerRow = ["guideId", "guide", "longSeq"]
    headerRow.extend(scoreNames)
    rows.insert(0, headerRow)
    return rows

def writeRow(ofh, row):
    " write list to file as tab-sep row "
    row = [str(x) for x in row]
    ofh.write("\t".join(row))
    ofh.write("\n")

def createBatchEffScoreTable(batchId):
    """ annotate all potential guides with efficiency scores and write to file.
    tab-sep file for easier debugging, no pickling
    """
    outFname = join(batchDir, batchId+".effScores.tab")
    seq, org, pam, position, extSeq = readBatchParams(batchId)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    guideRows = calcGuideEffScores(seq, extSeq, pam)
    guideFh = open(outFname, "w")
    for row in guideRows:
        writeRow(guideFh, row)
    guideFh.close()
    logging.info("Wrote eff scores to %s" % guideFh.name)

def readEffScores(batchId):
    " parse eff scores from tab sep file and return as dict pamId -> dict of scoreName -> value "
    effScoreFname = join(batchDir, batchId)+".effScores.tab"
    # old batches during transition time don't have this file yet, so make one now
    if not isfile(effScoreFname):
        createBatchEffScoreTable(batchId)

    seqToScores = {}
    for row in lineFileNext(open(effScoreFname)):
        scoreDict = {}
        rowDict = row._asdict()
        # the first three fields are the pamId, shortSeq, longSeq, they are not scores
        allScoreNames = row._fields[3:]
        for scoreName in allScoreNames:
            score = rowDict[scoreName]
            if "." in score:
                score = float(score)
            else:
                score = int(score)
            scoreDict[scoreName] = score
        seqToScores[row.guideId] = scoreDict
    return seqToScores



def makeVariants(seq):
    " generate all possible variants of sequence at 1bp-distance"
    seqs = []
    for i in range(0, len(seq)):
        for l in "ACTG":
            if l==seq[i]:
                continue
            newSeq = seq[:i]+l+seq[i+1:]
            seqs.append((i, seq[i], l, newSeq))
    return seqs

def expandIupac(seq):
    """ expand all IUPAC characters to nucleotides, returns list. 
    >>> expandIupac("NY")
    ['GC', 'GT', 'AC', 'AT', 'TC', 'TT', 'CC', 'CT']
    """
    # http://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
    d = {'A': 'A', 'C': 'C', 'B': 'CGT', 'D': 'AGT', 'G': 'G', \
        'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'GATC', 'S': 'CG', \
        'R': 'AG', 'T': 'T', 'W': 'AT', 'V': 'ACG', 'Y': 'CT', 'X': 'GATC'}
    seqs = []
    for i in product(*[d[j] for j in seq]):
       seqs.append("".join(i))
    return seqs

def writeBowtieSequences(inFaFname, outFname, pamPat):
    """ write the sequence and one-bp-distant-sequences + all possible PAM sequences to outFname 
    Return dict querySeqId -> querySeq and a list of all
    possible PAMs, as nucleotide sequences (not IUPAC-patterns)
    """
    ofh = open(outFname, "w")
    outCount = 0
    inCount = 0
    guideSeqs = {} # 20mer guide sequences
    qSeqs = {} # 23mer query sequences for bowtie, produced by expanding guide sequences
    allPamSeqs = expandIupac(pamPat)
    for seqId, seq in parseFastaAsList(open(inFaFname)):
        inCount += 1
        guideSeqs[seqId] = seq
        for pamSeq in allPamSeqs:
            # the input sequence + the PAM
            newSeqId = "%s.%s" % (seqId, pamSeq)
            newFullSeq = seq+pamSeq
            ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
            qSeqs[newSeqId] = newFullSeq

            # all one-bp mutations of the input sequence + the PAM
            for nPos, fromNucl, toNucl, newSeq in makeVariants(seq):
                newSeqId = "%s.%s.%d:%s>%s" % (seqId, pamSeq, nPos, fromNucl, toNucl)
                newFullSeq = newSeq+pamSeq
                ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
                qSeqs[newSeqId] = newFullSeq
                outCount += 1
    ofh.close()
    logging.debug("Wrote %d variants+expandedPam of %d sequences to %s" % (outCount, inCount, outFname))
    return guideSeqs, qSeqs, allPamSeqs

def applyModifStr(seq, modifStrs, strand):
    """ bowtie: given a list of pos:toNucl>fromNucl and a seq, return the original seq.
    position is 0-based
    >>> applyModifStr("ACAATAAGACATAAACATATCGG", "14:T>A,21:A>G,22:C>G".split(","), "+")
    'ACAATAAGACATAATCATATCAC'
    """
    seq = list(seq)
    for modifStr in modifStrs:
        #logging.debug( modifStr)
        pos, toFromNucl = modifStr.split(":")
        fromNucl, toNucl = toFromNucl.split(">")
        pos = int(pos)
        if strand=="-":
            fromNucl = revComp(fromNucl)
        seq[pos] = fromNucl
    return "".join(seq)
   
def parseRefout(tmpDir, guideSeqs, pamLen):
    """ parse all .map file in tmpDir and return as list of chrom,start,end,strand,guideSeq,tSeq
    """
    fnames = glob.glob(join(tmpDir, "*.map"))

    # while parsing, make sure we keep only the hit with the lowest number of mismatches
    # to the guide. Saves time when parsing.
    posToHit = {}
    hitBestMismCount = {}
    for fname in fnames:
        for line in open(fname):
           # s20+.17:A>G     -       chr8    26869044        CCAGCACGTGCAAGGCCGGCTTC IIIIIIIIIIIIIIIIIIIIIII 7       4:C>G,13:T>G,15:C>G
           guideIdWithMod, strand, chrom, start, tSeq, weird, someScore, alnModifStr = \
               line.rstrip("\n").split("\t")

           guideId = guideIdWithMod.split(".")[0]
           modifParts = alnModifStr.split(",")
           if modifParts==['']:
               modifParts = []
           mismCount = len(modifParts)
           hitId = (guideId, chrom, start, strand)
           oldMismCount = hitBestMismCount.get(hitId, 9999)
           if mismCount < oldMismCount:
               hit = (mismCount, guideIdWithMod, strand, chrom, start, tSeq, modifParts)
               posToHit[hitId] = hit

    ret = []
    for guideId, hit in posToHit.iteritems():
           mismCount, guideIdWithMod, strand, chrom, start, tSeq, modifParts = hit
           if strand=="-":
               tSeq = revComp(tSeq)
           guideId = guideIdWithMod.split(".")[0]
           guideSeq = guideSeqs[guideId]
           genomeSeq = applyModifStr(tSeq, modifParts, strand)
           start = int(start)
           bedRow = (guideId, chrom, start, start+GUIDELEN+pamLen, strand, guideSeq, genomeSeq) 
           ret.append( bedRow )

    return ret

def getEditDist(str1, str2):
    """ return edit distance between two strings of equal length 
    >>> getEditDist("HIHI", "HAHA")
    2
    """
    assert(len(str1)==len(str2))
    str1 = str1.upper()
    str2 = str2.upper()

    editDist = 0
    for c1, c2 in zip(str1, str2):
        if c1!=c2:
            editDist +=1
    return editDist

def processSubmission(faFname, genome, pam, bedFname, batchBase, batchId, queue):
    """ search fasta file against genome, filter for pam matches and write to bedFName 
    optionally write status updates to work queue.
    """
    doEffScoring = False
    if doEffScoring and not cpf1Mode:
        queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
        createBatchEffScoreTable(batchId)

    if genome=="noGenome":
        # skip off-target search
        if cpf1Mode:
            errAbort("Sorry, no efficiency score has been published yet for Cpf1.")
        open(bedFname, "w") # create a 0-byte file to signal job completion
        queue.startStep(batchId, "done", "Job completed")
        return

    if useBowtie:
        get_offtargets.findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pam, bedFname)
    else:
        get_offtargets.findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname)

    return bedFname

def lineFileNext(fh):
    """
        parses tab-sep file with headers as field names
        yields collection.namedtuples
        strips "#"-prefix from header line
    """
    line1 = fh.readline()
    line1 = line1.strip("\n").strip("#")
    headers = line1.split("\t")
    Record = namedtuple('tsvRec', headers)
   
    for line in fh:
        line = line.rstrip("\n")
        fields = line.split("\t")
        try:
            rec = Record(*fields)
        except Exception, msg:
            logging.error("Exception occured while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % repr(line))
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            #raise Exception("wrong field count in line %s" % line)
            continue
        # convert fields to correct data type
        yield rec

allGenomes = None

def readGenomes():
    " return list of all genomes supported "
    # global allGenomes
    if allGenomes:
        return allGenomes
    genomes = {}

    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")

    inFnames = []
    globalFname = join(genomesDir, "genomeInfo.all.tab")
    if isfile(globalFname):
        inFnames = [globalFname]
    else:
        for subDir in os.listdir(genomesDir):
            infoFname = join(genomesDir, subDir, "genomeInfo.tab")
            if isfile(infoFname):
                inFnames.append(infoFname)
    
    for infoFname in inFnames:
        for row in lineFileNext(open(infoFname)):
            # add a note to identify UCSC genomes
            if row.server.startswith("ucsc"):
                addStr="UCSC "
            else:
                addStr = ""
            genomes[row.name] = row.scientificName+" - "+row.genome+" - "+addStr+row.description

    genomes = genomes.items()
    genomes.sort(key=operator.itemgetter(1))
    allGenomes = genomes
    return allGenomes

def readBatchParams(batchId):
    """ given a batchId, return the genome, the pam, the input sequence and the
    chrom pos and extSeq, a 100bp-extended version of the input sequence.
    Returns None for pos if not found. """

    batchBase = join(batchDir, batchId)
    jsonFname = batchBase+".json"
    if isfile(jsonFname):
        params = json.load(open(jsonFname))
        # global batchName
        batchName = params["batchName"]
        return params["seq"], params["org"], params["pam"], params["posStr"], params["extSeq"]

    # FROM HERE UP TO END OF FUNCTION: legacy cold for old batches
    # remove in 2017
    inputFaFname = batchBase+".input.fa"
    if not isfile(inputFaFname):
        errAbort('Could not find the batch %s. We cannot keep Crispor runs for more than '
                'a few months. Please resubmit your input sequence via'
            ' <a href="crispor.py">the query input form</a>' % batchId)

    ifh = open(inputFaFname)
    ifhFields = ifh.readline().replace(">","").strip().split()
    if len(ifhFields)==2:
        genome, pamSeq = ifhFields
        position = None
    else:
        genome, pamSeq, position = ifhFields

    inSeq = ifh.readline().strip()

    ifh.seek(0)
    seqs = parseFasta(ifh)
    ifh.close()

    extSeq = None
    if "extSeq" in seqs:
        extSeq = seqs["extSeq"]

    # older batch files don't include a position yet
    if position==None:
        position = coordsToPosStr(*findBestMatch(genome, inSeq))

    return inSeq, genome, pamSeq, position, extSeq

def findAllPams(seq, pam):
    """ find all matches for PAM and return as dict startPos -> strand and a set
    of end positions
    """
    seq = seq.upper()
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet)
    return startDict, endSet

def newBatch(batchName, seq, org, pam, genomesDir,batchDir,skipAlign=False):
    """ obtain a batch ID and write seq/org/pam to their files.
    Return batchId, position string and a 100bp-extended seq, if possible.
    """
    batchId = makeTempBase(seq, org, pam, batchName)
    if skipAlign:
        chrom, start, end, strand = None, None, None, None
    else:
        chrom, start, end, strand = findBestMatch(org, seq,genomesDir)
    # define temp file names
    batchBase = join(batchDir, batchId)
    # save input seq, pamSeq, genome, position for primer design later
    batchJsonName = batchBase+".json"
    posStr = coordsToPosStr(chrom, start, end, strand)

    ofh = open(batchJsonName, "w")
    #ofh.write(">%s %s %s\n%s\n" % (org, pam, posStr, batchName, seq))
    batchData = {}
    batchData["org"] = org
    batchData["pam"] = pam
    batchData["posStr"] = posStr
    batchData["batchName"] = batchName
    batchData["seq"] = seq

    # try to get a 100bp-extended version of the input seq
    extSeq = None
    if chrom!=None:
        extSeq = extendAndGetSeq(org, genomesDir,chrom, start, end, strand)
        #if extSeq!=None:
            #ofh.write(">extSeq\n%s\n" % (extSeq))
    batchData["extSeq"] = extSeq

    json.dump(batchData, ofh)
    ofh.close()
    return batchId, posStr, extSeq

def readDbInfo(org):
    " return a dbInfo object with the columsn in the genomeInfo.tab file "
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")
    if not isfile(infoFname):
        return None
    dbInfo = lineFileNext(open(infoFname)).next()
    return dbInfo


def getOfftargets(seq, org, pam, batchId,batchDir, startDict, queue):
    """ write guides to fasta and run bwa or use cached results.
    Return name of the BED file with the matches.
    Write progress status updates to queue object.
    """
    batchBase = join(batchDir, batchId)
    otBedFname = batchBase+".bed"
    # write potential PAM sites to file 
    faFname = batchBase+".fa"
    writePamFlank(seq, startDict, pam, faFname)
    processSubmission(faFname, org, pam, otBedFname, batchBase, batchId, queue)

    return otBedFname

def getSeq(db, posStr):
    """
    given a database name and a string with the position as chrom:start-end, return the sequence as
    a string.
    """
    chrom, start, end, strand =  parsePos(posStr)
    if end-start > MAXSEQLEN and db!="noGenome":
        errAbort("Input sequence range too long. Please retry with a sequence range shorter than %d bp." % MAXSEQLEN)
    genomeDir = genomesDir # pull in global var
    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    binPath = join(binDir, "twoBitToFa")
    cmd = [binPath, twoBitFname, "-seq="+chrom, "-start="+str(start), "-end="+str(end), "stdout"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()
    # remove fasta header line
    lines = seqStr.splitlines()
    if len(lines)>0:
        lines.pop(0)
    seq = "".join(lines)
    if len(seq) < 23:
        errAbort("Sorry, the sequence range %s on genome %s is not longer than 23bp. To find a valid CRISPR/Cas9 site, one needs at least a 23bp long sequence." % (db, posStr))
    return seq


def readVarDbs(db):
    """ find all possible variant VCFs and return as list of (shortLabel, fname, label, hasAF) 
    hasAF = file has the AF field (allele frequency). Means that the UI
    will show the "frequency filter" button.
    """
    # parse the descriptions of the VCF files
    # descriptions are optional
    labelFname = join(genomesDir, db, "vcfDescs.txt")
    ret = []
    if isfile(labelFname):
        for line in open(labelFname):
            if line.startswith("#"):
                continue
            fields = string.split(line.rstrip("\n"), "\t")
            if len(fields)==4:
                shortLabel, fname, desc, hasAF = fields
            else:
                errAbort("not four fields in vcfDescs.txt: %s" % fields)

            fpath = join(genomesDir, db, fname)
            if not isfile(fpath):
                print "Error: Cannot find VCF file %s" % fpath
                continue
            hasAF = (hasAF=="1")
            ret.append( (shortLabel, fname, desc, hasAF) )
    return ret
        
def parseVcfInfo(info):
    " parse a VCF info string and return as a dict key=val "
    fs = info.split(";")
    ret = {}
    for field in fs:
        parts = string.split(field, "=", maxsplit=1)
        if len(parts)==2:
            key, val = parts
        elif len(parts)==1:
            key = parts[0]
            val = True
        ret[key] = val
    return ret

def findVariantsInRange(vcfFname, chrom, start, end, strand, minFreq):
    """ find variants that overlap the position.
    varDb is a tuple of label, vcfFname
    return as a dict relative position -> (chrom, pos, refAllele, altAllele, list of info-dicts)
    special position is "label" which is the label of the variant db
    """
    minFreq = float(minFreq)
    seqLen = end-start
    if not isfile(vcfFname):
        errAbort("%s not found" % vcfFname)
    tb = tabix.open(vcfFname)
    chrom = chrom.replace("chr","")
    records = tb.query(chrom, start+1, end) # VCF is 1-based

    varDict = defaultdict(list)
    for rec in records:
        chrom, varPos, varId, refAll, altAllStrList, qual, filterFlag, info = rec[:8]
        infoDict = parseVcfInfo(info)
        altAllList = altAllStrList.split(",")
        if "AF" in infoDict:
            afList = infoDict["AF"].split(",")
        else:
            afList = [None] * len(altAllList)
            
        for altAll, allFreq in zip(altAllList, afList):
            if minFreq is not None and allFreq is not None:
                allFreq = float(allFreq)
                if not allFreq > minFreq:
                    continue
                
            attribs = {}
            if allFreq!=None:
                attribs["freq"] = allFreq
            relPos = int(varPos)-1-start
            if strand=="-":
                relPos = seqLen - relPos - len(refAll)
                refAll = revComp(refAll)
                altAlls = []
                for altAll in altAll.split(","):
                    altAlls.append(revComp(altAll))
                altAll = ",".join(altAlls)
            if varId != ".":
                attribs["varId"] = varId
            varInfo = (chrom, varPos, refAll, altAll, attribs)
            varDict[relPos].append(varInfo)

    return varDict

def crisprSearch(params):
    " do crispr off target search and eff. scoring "
    # retrieve sequence if not provided
    if "pos" in params and not "seq" in params:
        params["seq"] = getSeq(params["org"], params["pos"])

    if "batchId" in params:
        # if we're getting only the batchId, extract the parameters from the batch
        # this allows a stable link to a batch that is done
        batchId = params["batchId"]
        seq, org, pam, position, extSeq = readBatchParams(batchId)
        seq, warnMsg = cleanSeq(seq, org)
    else:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        newBatchName = params.get("name", "")

        # the "seq" parameter can contain a chrom:start-end position instead of the sequence.
        if re.match(" *[a-zA-Z0-9_-]+: *[0-9, ]+ *- *[0-9,]+ *", seq):
            seq = getSeq(params["org"], seq.replace(" ","").replace(",",""))

        seq, warnMsg = cleanSeq(seq, org)
        batchId, position, extSeq = newBatch(newBatchName, seq, org, pam) 
    GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(pam)

    # check if minFreq was specified
    minFreq = params.get("minFreq", "0.0")
    try:
        minFreq = float(minFreq)
    except ValueError:
        errAbort("minFreq has to be a floating point number")

    varDb = params.get("varDb", None)

    if len(warnMsg)!=0:
        print warnMsg+"<p>"

    batchBase = join(batchDir, batchId)

    # read genome info tab file into memory
    dbInfo = readDbInfo(org)

    # search PAMs
    uppSeq = seq.upper()
    startDict, endSet = findAllPams(uppSeq, pam)
    otBedFname = getOfftargets(uppSeq, org, pam, batchId, startDict, None)
    if position=='?':
        printQueryNotFoundNote(dbInfo)
    else:
        genomePosStr = ":".join(position.split(":")[:2])
        chrom, start, end, strand = parsePos(position)
        start = str(int(start)+1)
        oneBasedPosition = "%s:%s-%s" % (chrom, start, end)

        print "<div class='title'><em>"
        if batchName!="":
            print batchName+":"

        print "%s (%s)</em>, " % (dbInfo.scientificName, dbInfo.name)
        print '<span style="text-decoration:underline">'
        #mouseOver = "link to UCSC,Ensembl or Gbrowse Genome Browser"
        mouseOver = None
        if dbInfo.server=="manual":
            mouseOver = "no genome browser link available for this organism"
        if strand=="+":
            print " forward genomic strand"
        else:
            print " reverse genomic strand"
        print "</div>"
        #print " (link to Genome Browser)</div>"

    otMatches = parseOfftargets(otBedFname)
    effScores = readEffScores(batchId)
    sortBy = (params.get("sortBy", None))
    guideData, guideScores, hasNotFound = mergeGuideInfo(uppSeq, startDict, pam, otMatches, position, effScores, sortBy)

    if len(guideScores)==0:
        print "Found no possible guide sequence. Make sure that your input sequence is long enough and contains at least one match to the PAM motif %s." % pam
        print '<br><a class="neutral" href="crispor.py">'
        print '<div class="button" style="margin-left:auto;margin-right:auto;width:150;">New Query</div></a>'
        return

    if hasNotFound and not position=="?":
        print('<div style="text-align:left"><strong>Note:</strong> At least one of the possible guide sequences was not found in the genome. ')
        print("If you pasted a cDNA sequence, note that sequences with score 0, e.g. splice junctions, are not in the genome, only in the cDNA and are not usable as CRISPR guides.</div><br>")

    chrom, start, end, strand = parsePos(position)

    # get list of variant databases
    varLabel = None
    varDbs = readVarDbs(org)

    if len(varDbs)>0:
        if varDb is None:
            varDb = varDbs[0][1]

        # pull out label of the variant database
        varLabel = None
        for shortLabel, varKey, lab, hasAF in varDbs:
            if varKey==varDb:
                varLabel = lab
                break
        if varLabel is None:
            errAbort("variant DB %s was not found in vcfDescs.txt" % varDb)

        vcfFname = join(genomesDir, org, varDb)
        varDict = findVariantsInRange(vcfFname, chrom, start, end, strand, minFreq)
        varDict["label"] = varLabel
    else:
        varDict = None
        varLabel = None

    
pamIdRe = re.compile(r's([0-9]+)([+-])g?([0-9]*)')

def intToExtPamId(pamId):
    " convert the internal pam Id like s20+ to the external one, like 21Forw "
    pamPos, strand, rest = pamIdRe.match(pamId).groups()
    if strand=="+":
        strDesc = 'forw'
    else:
        strDesc = 'rev'
    guideDesc = str(int(pamPos)+1)+strDesc
    return guideDesc

def concatGuideAndPam(guideSeq, pamSeq):
    " return guide+pam or pam+guide, depending on cpf1Mode "
    if cpf1Mode:
        return pamSeq+guideSeq
    else:
        return guideSeq+pamSeq

def iterGuideRows(guideData, addHeaders=False):
    "yield rows from guide data. Need to know if for Cpf1 or not "
    headers = list(tuple(guideHeaders)) # make a copy of the list
    for scoreName in scoreNames:
        headers.append(scoreDescs[scoreName][0]+"EffScore")
    if addHeaders:
        yield headers

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow

        otCount = 0
        if otData!=None:
            otCount = len(otData)

        guideDesc = intToExtPamId(pamId)

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)
        row = [guideDesc, fullSeq, guideScore, guideCfdScore, otCount, ontargetDesc]
        for scoreName in scoreNames:
            row.append(effScores.get(scoreName, "NotEnoughFlankSeq"))
        row = [str(x) for x in row]
        yield row

def iterOfftargetRows(guideData, addHeaders=False):
    " yield bulk offtarget rows for the tab-sep download file "
    if addHeaders:
        headers = list(offtargetHeaders) # clone list
        yield headers

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, \
            ontargetDesc, subOptMatchCount = guideRow

        otCount = 0

        if otData!=None:
            otCount = len(otData)
            for otSeq, mitScore, cfdScore, editDist, pos, gene, alnHtml in otData:
                gene = gene.replace(",", "_").replace(";","-")
                chrom, start, end, strand = parsePos(pos)
                guideDesc = intToExtPamId(pamId)
                mismStr = highlightMismatches(guideSeq, otSeq, len(pamSeq))
                fullSeq = concatGuideAndPam(guideSeq, pamSeq)
                row = [guideDesc, fullSeq, otSeq, mismStr, editDist, mitScore, cfdScore, chrom, start, end, strand, gene]
                row = [str(x) for x in row]
                yield row

def xlsWrite(rows, title, outFile, colWidths, fileFormat):
    """ given rows, writes a XLS binary stream to outFile, if xlwt is available 
    Otherwise writes a tab-sep file.
    colWidths is a list of widths of columns, in Arial characters.
    """
    if xlwtLoaded and not fileFormat=="tsv":
        charSize = 269 # see http://reliablybroken.com/b/2011/10/widths-heights-with-xlwt-python/
        wb = xlwt.Workbook()
        ws = wb.add_sheet(title)

        for rowCount, row in enumerate(rows):
            for colCount, col in enumerate(row):
                if col.isdigit():
                    col = int(col)
                ws.write(rowCount, colCount, col)

        # set sizes in characters per column
        for colId, colWidth in enumerate(colWidths):
            ws.col(colId).width = charSize*colWidth

        wb.save(outFile)
    else:
        for row in rows:
            outFile.write("\t".join(row))
            outFile.write("\n")
    outFile.flush()

def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs.append( (seqId, "".join(parts)) )
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs.append( (seqId, "".join(parts)) )
    return seqs

def parseFasta(fileObj):
    " parse a fasta file, where each seq is on a single line, return dict id -> seq "
    seqs = {}
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs[seqId]  = "".join(parts)
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs[seqId]  = "".join(parts)
    return seqs

def coordsToPosStr(chrom, start, end, strand):
    " convert coords to a string "
    if chrom==None:
        return "?"
    locStr = "%s:%d-%d:%s" % (chrom, start, end, strand)
    return locStr

def findBestMatch(genome, seq,genomesDir):
    """ find best match for input sequence from batchId in genome and return as
    a string "chrom:start-end:strand or None if not found "
    """
    if genome=="noGenome":
        return None, None, None, None

    # write seq to tmp file
    tmpFaFh = tempfile.NamedTemporaryFile(prefix="crisporBestMatch", suffix=".fa")
    tmpFaFh.write(">tmp\n%s" % seq)
    tmpFaFh.flush()
    logging.debug("seq: %s" % open(tmpFaFh.name).read())
    faFname = tmpFaFh.name

    # create temp SAM file
    tmpSamFh = tempfile.NamedTemporaryFile(prefix="crisporBestMatch", suffix=".sam")
    samFname = tmpSamFh.name

    genomeDir = genomesDir # make local var
    cmd = "$BIN/bwa bwasw -T 20 %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(samFname)s" % locals()
    get_offtargets.runCmd(cmd)

    chrom, start, end = None, None, None
    for l in open(samFname):
        if l.startswith("@"):
            continue
        logging.debug("%s" % l)
        l = l.rstrip("\n")
        fs = l.split("\t")
        qName, flag, rName, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fs[:11]
        if (int(flag) and 2) == 2:
            strand = "-"
        else:
            strand = "+"
        if not re.compile("[0-9]*").match(cigar):
            continue
        if cigar=="*":
            logging.debug("No best match found")
            return None, None, None, None
            #errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        # XX why do we get soft clipped sequences from BWA? repeats?
        cleanCigar = cigar.replace("M","").replace("S", "")
        if not cleanCigar.isdigit():
            logging.debug("Best match found, but cigar string was %s" % cigar)
            return None, None, None, None
        matchLen = int(cleanCigar)
        chrom, start, end =  rName, int(pos)-1, int(pos)-1+matchLen # SAM is 1-based
        #print chrom, start, end, strand

    # delete the temp files
    tmpSamFh.close()
    tmpFaFh.close()
    logging.debug("Found best match at %s:%d-%d:%s" % (chrom, start, end, strand))
    return chrom, start, end, strand

def getGenomeSeqs(genome, coordList):
    """ return dict of genome sequences,
    coordList has format (chrom, start, end, name)
    returns list (chrom, start, end, name, seq)
    """
    binFname = join(binDir, "twoBitToFa")
    twoBitPath = "genomes/%(genome)s/%(genome)s.2bit" % locals()
    twoBitPath = abspath(twoBitPath)

    tbf = twobitreader.TwoBitFile(twoBitPath)
    seqs = []
    for coordTuple in coordList:
        chrom, start, end, name = coordTuple
        seqs.append((chrom, start, end, name, tbf[chrom][start:end]) )
    return seqs


def mergeParamDicts(params, changeParams):
    """ changeParams is a dict that can override elements in params.
    if value==None in changeParams, the whole element will get removed.
    if onlyParams is set, only copy over the keys in onlyParams (a list)
    """
    newParams = {}
    newParams.update(params)
    newParams.update(changeParams)
    for key, val in changeParams.iteritems():
        if val==None:
            del newParams[key]
    return newParams




def findGuideSeq(inSeq, pam, pamId):
    """ given the input sequence and the pamId, return the guide sequence,
    the sequence with the pam and its strand.
    """
    startDict, endSet = findAllPams(inSeq, pam)
    pamInfo = list(flankSeqIter(inSeq, startDict, pam, False))
    for guidePamId, pamStart, guideStart, guideStrand, guideSeq, pamSeq in pamInfo:
        if guidePamId!=pamId:
            continue

        guideSeqWPam = guideSeq+pamSeq
        # prettify guideSeqWPam to highlight the PAM
        if cpf1Mode:
            guideSeqHtml = "<i>%s</i> %s" % \
                (guideSeqWPam[:len(pam)].upper(), guideSeqWPam[len(pam):].upper())
        else:
            guideSeqHtml = "%s <i>%s</i>" % \
                (guideSeqWPam[:-len(pam)].upper(), guideSeqWPam[-len(pam):].upper())

        guideEnd = guideStart + GUIDELEN
        return guideSeq, pamSeq, guideSeqWPam , guideStrand, guideSeqHtml, \
                guideStart, guideEnd
    errAbort("pamId %s not found? This is a bug." % pamId)

def findOntargetPos(batchBase, pamId, position):
    " find position of guide sequence in genome at MM0 "
    otBedFname = batchBase+".bed"
    otMatches = parseOfftargets(otBedFname)
    if pamId not in otMatches or 0 not in otMatches[pamId]:
        errAbort("No perfect match found for guide sequence in the genome. Cannot design primers for a non-matching guide sequence.<p>Are you sure you have selected the right genome? <p> If you have selected the right genome and entered a cDNA as the query sequence, please note that sequences that overlap a splice site are not part of the genome and cannot be used as guide sequences.")

    matchList = otMatches[pamId][0] # = get all matches with 0 mismatches
    if len(matchList)>1:
        targetChrom, targetStart, targetEnd, strand = parsePos(position)

        filtMatch = None
        # search for off-target that is the on-target
        for match in matchList:
            # example match: ('scaffold_1', 578, 601, 'TATTGGATTGGTCCAATCGTTGG', '-', 'ex', 'GSADVT00000001001', 293) 
            chrom, start, end = match[:3]
            if chrom==targetChrom and start>=targetStart and end<=targetEnd:
                filtMatch = match
                break

        if filtMatch is None:
            errAbort("Multiple matches for this guide, but no single match is within the target sequence? Please contact us at %s, this looks like a bug or at the very least an issue that was not tested." % contactEmail)

        chrom, start, end = filtMatch[:3]
        gene = filtMatch[6]
        print("<strong>Warning</strong>: Found multiple perfect matches for this guide sequence in the genome. For the PCR, we are using the on-target match in the input sequence %s:%d-%d (gene: %s), but this guide will not be specific. Is this a polyploid organism? Try selecting another guide sequence or email %s to discuss your strategy or modifications to this software.<p>" % (chrom, start+1, end, gene, contactEmail))

        matchList = [filtMatch]

    # global batchName
    batchName = batchName.replace(" ", "_")

    chrom, start, end, seq, strand, segType, segName, x1Count = \
        matchList[0]
    start = int(start)
    end = int(end)
    return chrom, start, end, strand


def runTests():
    guideSeq = "CTCTTTACGCAGAGGGATGT"
    testRes = {"ATTTTTATGCAGAGTGATGT":     0.4, 
               "TTCTTTACCCGGAGGGATGA": 0.2, 
               "CTGTTTACACACAGGGATTT": 0.2, 
               "CTCTCTGTGCAGATGGATGT": 0.1, 
               "ATCTTAAAGCAGATGGATGT": 0.1, 
               "CTCTTTCCGCAGAGGCTTGT": 0.1, 
               "CTCGTAGCGCAGAGGGAGGT": 0.1, 
               "CTCTTTAAAGAGATGGATGT": 0.1, 
               "CACTTCACTCAGAGGCATGT": 0.1, 
               "CTTTTTTCTCAGAAGGATGT": 0.1, 
               "CTCTTTACACAGAGAGACGT": 0.1, 
               "CTCTTTTCTCAGAGAGATGG": 0.1, 
               "CTATTTACCCAAATGGATGT": 0.1, 
               "CTCTTTGCACAGGGGGAAGT": 0, 
               "CTCTTTGCACAGGGGGAAGT": 0, 
               "CTCTTCACACAGAGGAATGA": 0, 
               "CTCTTTCCACAGGGGAATGT": 0 }

    testRes2 = {
       "GAGTCTAAGCAGAAGAAGAA":     2.2,
       "GAGTCCTAGCAGGAGAAGAA": 1.8,
       "GAGAGCAAGCAGAAGAAGAA": 1.6,
       "GAGTACTAGAAGAAGAAGAA": 1.6,
       "ACGTCTGAGCAGAAGAAGAA": 1.5,
       "GCGACAGAGCAGAAGAAGAA": 1.5,
       "GAGTAGGAGGAGAAGAAGAA": 1.4,
       "GATGCCGTGAAGAAGAAGAA": 1.3,
       "GATTCCTACCAGAAGAAGAA": 1,
       "GAATCCAAGCAGAAGAAGAG": 1,
       "AAGTACTGGCAGAAGAAGAA": 0.9,
       "AGGTGCTAGCAGAAGAAGAA": 0.9,
       "GGGGCCAGGCAGAAGAAGAA": 0.9,
       "ATGTGCAAGCAGAAGAAGAA": 0.9,
       "ACCTCCCAGCAGAAGAAGAA": 0.9,
       "CCCTCCCAGCAGAAGAAGAA": 0.9,
       "TCATCCAAGCAGAAGAAGAA": 0.9,
       "TTCTCCAAGCAGAAGAAGAA": 0.9,
       "GGTGCCAAGCAGAAGAAGAA": 0.9,
       "GCACCCCAGCAGAAGAAGAA": 0.9,
       "CAGTCCAGGAAGAAGAAGAA": 0.9,
       "AAGCCCAAGGAGAAGAAGAA": 0.9,
       "CACTCCAAGTAGAAGAAGAA": 0.9,
       "GAGTCCGGGAAGGAGAAGAA": 0.9,
       "GGTTCCCAGGAGAAGAAGAA": 0.9,
       "AAGTCTGAGCACAAGAAGAA": 0.9,
       "GAGGACAAGAAGAAGAAGAA": 0.9,
       "GTCTGCGATCAGAAGAAGAA": 0.8,
       "GGTTCTGTGCAGAAGAAGAA": 0.8,
       "AGGTGGGAGCAGAAGAAGAA": 0.8,
       "AAGAGCGAGCGGAAGAAGAA": 0.8,
       "CAATTTGAGCAGAAGAAGAA": 0.8,
       "AATACAGAGCAGAAGAAGAA": 0.8,
       "CAAACGGAGCAGAAGAAGAA": 0.8,
       "AAGTGAGAGTAGAAGAAGAA": 0.8,
       "AAGTAGGAGAAGAAGAAGAA": 0.8,
       "AAGTTGGAGAAGAAGAAGAA": 0.8,
       "CAGGCTGAGAAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "AAGTGGGAGGAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "TCTTCCGAGCTGAAGAAGAA": 0.8,
       "GCGGCCGATGAGAAGAAGAA": 0.8,
       "GCGTCCGCCAAGAAGAAGAA": 0.8,
       "GCTCCTGAGCAGAAGAAGAA": 0.8,
       "CACTCTGAGGAGAAGAAGAA": 0.8,
       "GTGTGGGAGGAGAAGAAGAA": 0.8,
       "GGGTAAGAGTAGAAGAAGAA": 0.8
    }
    #for seq, expScore in testRes.iteritems():
        #score = calcHitScore(guideSeq, seq)
        #print score, "%0.1f" % score, expScore

    guideSeq = "GAGTCCGAGCAGAAGAAGAA"
    for seq, expScore in testRes2.iteritems():
        score = calcHitScore(guideSeq, seq)
        #print score, "%0.1f" % score, expScore
    


def delBatchDir():
    " called at program exit, for command line mode "
    delTmpDirs() # first remove any subdirs
    if not isdir(batchDir):
        return
    logging.debug("Deleting dir %s" % batchDir)
    fnames = glob.glob(join(batchDir, "*"))
    if len(fnames)>50:
        raise Exception("cowardly refusing to remove many temp files")
    for fname in fnames:
        os.remove(fname)
    os.removedirs(batchDir)

tmpDirsDelExit = []

def delTmpDirs():
    " signal handler at program exit, to remove registered tmp dirs "
    # global tmpDirsDelExit
    logging.debug("Removing tmpDirs: %s" % ",".join(tmpDirsDelExit))
    for tmpDir in tmpDirsDelExit:
        if isdir(tmpDir):
            shutil.rmtree(tmpDir)
    tmpDirsDelExit = []



def handleOptions(options):
    " set glpbal vars based on options "
    from constants import DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,useBowtie
    if options.test:
        runTests()
        import doctest
        doctest.testmod()
        sys.exit(0)

    if options.debug:
        DEBUG = True

    if options.noEffScores:
        doEffScoring = False

    # handle the alignment/filtering options
    if options.maxOcc != None:
        MAXOCC = options.maxOcc

    if options.minAltPamScore!=None:
        ALTPAMMINSCORE = options.minAltPamScore

    if options.mismatches:
        maxMMs = options.mismatches

    if options.bowtie:
        useBowtie = True

    if options.pam:
        GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(options.pam)
    argument_list = [DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,
                            useBowtie,GUIDELEN,cpf1Mode,addGenePlasmids]
    
    return argument_list

    
def mainCommandLine():
    " main entry if called from command line "
    # global commandLineMode
    commandLineMode = True

    args, options = parser_functions.parseArgs()
    argument_list = handleOptions(options)
    [DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,
    useBowtie,GUIDELEN,cpf1Mode,addGenePlasmids] = argument_list
    org, inSeqFname, outGuideFname = args

    skipAlign = False
    if options.skipAlign:
        skipAlign = True

    # different genomes directory?
    if options.genomeDir != None:
        genomesDir = options.genomeDir

    # get sequence
    seqs = parseFasta(open(inSeqFname))

    # make a directory for the temp files
    # and put it into a global variable, so all functions will use it
    # global batchDir
    batchDir = tempfile.mkdtemp(dir=TEMPDIR, prefix="crispor")
    logging.debug("Created directory %s" % batchDir)
    if options.debug:
        logging.info("debug-mode, temporary directory %s will not be deleted" % batchDir)
    else:
        atexit.register(delBatchDir)

    # prepare output files
    guideFh = open(join(batchDir, "guideInfo.tab"), "w")
    guideFh.write("\t".join(guideHeaders)+"\n")
    if options.offtargetFname:
        offtargetFh = open(join(batchDir, "offtargetInfo.tab"), "w")
        offtargetFh.write("\t".join(offtargetHeaders)+"\n")

    # putting multiple sequences into the input file is possible
    # but very inefficient. Rather separate them with a stretch of 10 Ns
    # as explained the docs
    for seqId, seq in seqs.iteritems():
        seq = seq.upper()
        logging.info("running on sequence ID '%s'" % seqId)
        # get the other parameters and write to a new batch
        seq = seq.upper()
        pam = options.pam
        batchId, position, extSeq = newBatch(seqId, seq, org, pam,genomesDir,batchDir,skipAlign)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        if position=="?":
            logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pam)
        otBedFname = getOfftargets(seq, org, pam, batchId, batchDir,startDict, ConsQueue())
        otMatches = parseOfftargets(otBedFname)

        if options.noEffScores or cpf1Mode:
            effScores = {}
        else:
            effScores = readEffScores(batchId)

        guideData, guideScores, hasNotFound = \
            mergeGuideInfo(seq, startDict, pam, otMatches, position, effScores)

        for row in iterGuideRows(guideData):
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if options.offtargetFname:
            for row in iterOfftargetRows(guideData):
                offtargetFh.write("\t".join(row))
                offtargetFh.write("\n")

    guideFh.close()
    shutil.move(guideFh.name, outGuideFname)
    logging.info("guide info written to %s" % outGuideFname)

    if options.offtargetFname:
        offtargetFh.close()
        shutil.move(offtargetFh.name, options.offtargetFname)
        logging.info("off-target info written to %s" % options.offtargetFname)

    if not options.debug:
       shutil.rmtree(batchDir)


def main():
    mainCommandLine()
main()
