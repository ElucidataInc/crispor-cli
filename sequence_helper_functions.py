from os.path import join
from constants import baseDir
from collections import defaultdict, namedtuple
    
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
    cpf1Mode = (pamSeq=='TTTN')
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

def concatGuideAndPam(guideSeq, pamSeq):
    " return guide+pam or pam+guide, depending on cpf1Mode "
    cpf1Mode= (pamSeq=='TTTN')
    if cpf1Mode:
        return pamSeq+guideSeq
    else:
        return guideSeq+pamSeq

