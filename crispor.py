#!/usr/bin/env python2.7
# the tefor crispr tool
# can be run as a CGI or from the command line

# OOF scores are WRONG for Cpf1! -> where is the cut site?

# python std library
import subprocess, tempfile, logging, atexit, glob, shutil
import Cookie, time, sys, cgi, re, random, platform, os
import hashlib, base64, string, operator, urllib, sqlite3, time
import traceback,  pwd, pickle
import common_functions
import parser_functions
import get_offtargets
import batch_functions
import sequence_helper_functions
import scoring_functions

from bio_functions import *
from fasta_functions import *
from common_functions import *
from basic_imports import *
from constants import *
from collections import defaultdict, namedtuple
from datetime import datetime
from itertools import product
from os.path import join, isfile,  dirname, isdir, abspath

 

concatGuideAndPam = sequence_helper_functions.concatGuideAndPam


# ====== FUNCTIONS =====
contentLineDone = False




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
                mitScore = scoring_functions.calcHitScore(guideNoPam, otSeqNoPam)
                # CFD score must include the PAM
                cfdScore = scoring_functions.calcCfdScore(guideSeq, otSeq)
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
            guideScore = scoring_functions.calcMitGuideScore(sum(mitOtScores))
            guideCfdScore = scoring_functions.calcMitGuideScore(sum(cfdScores))

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


def mergeGuideInfo(seq, startDict, pamPat, otMatches, inputPos, effScores, sortBy=None):
    """
    merges guide information from the sequence, the efficiency scores and the off-targets.
    creates rows with fields:


    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "effScore", "mhScore", "oofScore" or "pos"
    """
    allEnzymes = sequence_helper_functions.readEnzymes()

    guideData = []
    guideScores = {}
    hasNotFound = False

    pamSeqs = list(get_offtargets.flankSeqIter(seq, startDict, pamPat, True))

    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = concatGuideAndPam(guideSeq, pamSeq)
            mutEnzymes = sequence_helper_functions.matchRestrEnz(allEnzymes, guideSeq, pamSeq)
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




def debug(msg):
    logging.debug(msg)
    if DEBUG:
        print msg



class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))




allGenomes = None




def findAllPams(seq, pam):
    """ find all matches for PAM and return as dict startPos -> strand and a set
    of end positions
    """
    seq = seq.upper()
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, common_functions.revComp(pam), "-", startDict, endSet)
    return startDict, endSet




    
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
                mismStr = get_offtargets.highlightMismatches(guideSeq, otSeq, len(pamSeq))
                fullSeq = concatGuideAndPam(guideSeq, pamSeq)
                row = [guideDesc, fullSeq, otSeq, mismStr, editDist, mitScore, cfdScore, chrom, start, end, strand, gene]
                row = [str(x) for x in row]
                yield row

def mainCommandLine():
    " main entry if called from command line "
    args, options = parser_functions.parseArgs()
    argument_list = parser_functions.handleOptions(options)
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
    seqs = batch_functions.parseFasta(open(inSeqFname))

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
        batchId, position, extSeq = batch_functions.newBatch(seqId, seq, org, pam,genomesDir,batchDir,skipAlign)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        if position=="?":
            logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pam)
        process_parameters = [doEffScoring,useBowtie ,cpf1Mode]
        otBedFname = get_offtargets.getOfftargets(seq, org, pam, batchId, batchDir,startDict, ConsQueue(),process_parameters)
        otMatches = get_offtargets.parseOfftargets(otBedFname)

        if options.noEffScores or cpf1Mode:
            effScores = {}
        else:
            effScores = readEffScores(batchDir,batchId)

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
