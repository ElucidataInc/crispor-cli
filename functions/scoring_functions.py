import batch_functions
import common_functions
import crisporEffScores
import fasta_functions
import json
import logging
import pickle
import re
import sequence_helper_functions

from collections import defaultdict, namedtuple
from os.path import join,dirname,isfile




# --- START OF SCORING ROUTINES 

# MIT offtarget scoring

# aka Matrix "M"
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]


def get_mm_pam_scores():
    """
    """
    dataDir = join(dirname(__file__), '../CFD_Scoring')
    mm_scores = pickle.load(open(join(dataDir, 'mismatch_score.pkl'),'rb'))
    pam_scores = pickle.load(open(join(dataDir, 'pam_scores.pkl'),'rb'))
    return (mm_scores,pam_scores)

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)



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

    startDict, endSet = fasta_functions.findAllPams(seq, pam)
    pamInfo = list(fasta_functions.flankSeqIter(seq, startDict, pam, False))

    guideIds = []
    guides = []
    longSeqs = []
    for pamId, startPos, guideStart, strand, guideSeq, pamSeq in pamInfo:
        guides.append(guideSeq+pamSeq)
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam),cpf1Mode,GUIDELEN)
        longSeq = sequence_helper_functions.getExtSeq(seq, gStart, gEnd, strand, 50-GUIDELEN, 50, extSeq)
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

def createBatchEffScoreTable(batchDir,batchId):
    """ annotate all potential guides with efficiency scores and write to file.
    tab-sep file for easier debugging, no pickling
    """
    outFname = join(batchDir, batchId+".effScores.tab")
    seq, org, pam, position, extSeq = readBatchParams(batchDir,batchId)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    guideRows = calcGuideEffScores(seq, extSeq, pam)
    guideFh = open(outFname, "w")
    for row in guideRows:
        writeRow(guideFh, row)
    guideFh.close()
    logging.info("Wrote eff scores to %s" % guideFh.name)

def readEffScores(batchDir,batchId):
    " parse eff scores from tab sep file and return as dict pamId -> dict of scoreName -> value "
    effScoreFname = join(batchDir, batchId)+".effScores.tab"
    # old batches during transition time don't have this file yet, so make one now
    if not isfile(effScoreFname):
        createBatchEffScoreTable(batchDir,batchId)

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

def readBatchParams(batchDir,batchId):
    """ given a batchId, return the genome, the pam, the input sequence and the
    chrom pos and extSeq, a 100bp-extended version of the input sequence.
    Returns None for pos if not found. """

    batchBase = join(batchDir, batchId)
    jsonFname = batchBase+".json"
    if isfile(jsonFname):
        params = json.load(open(jsonFname))
        global batchName
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
        position = batch_functions.coordsToPosStr(*findBestMatch(genome, inSeq))

    return inSeq, genome, pamSeq, position, extSeq


