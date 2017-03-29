import common_functions

def findPat(seq, pat):
    """ yield positions where pat matches seq, stupid brute force search 
    """
    seq = seq.upper()
    pat = pat.upper()
    for i in range(0, len(seq)-len(pat)+1):
        #print "new pos", i, seq[i:i+len(pat)],"<br>"
        found = True
        for x in range(0, len(pat)):
            #print "new step", x, "<br>"
            if pat[x]=="N":
                #print "N","<br>"
                continue
            seqPos = i+x
            if seqPos == len(seq):
                found = False
                break
            if not matchNuc(pat[x], seq[seqPos]):
            #if not patMatch(seq[seqPos], pat[x]):
                #print i, x, pat[x], seq[seqPos], "no match<br>"
                found = False
                break
            #print "match", i, x, found, "<br>"
        if found:
            #print "yielding", i, "<br>"
            yield i

def rndSeq(seqLen):
    " return random seq "
    seq = []
    alf = "ACTG"
    for i in range(0, seqLen):
        seq.append(alf[random.randint(0,3)])
    return "".join(seq)

def cleanSeq(seq, db):
    """ remove fasta header, check seq for illegal chars and return (filtered
    seq, user message) special value "random" returns a random sequence.
    """
    #print repr(seq)
    if seq.startswith("random"):
        seq = rndSeq(800)
    lines = seq.strip().splitlines()
    #print "<br>"
    #print "before fasta cleaning", "|".join(lines)
    if len(lines)>0 and lines[0].startswith(">"):
        line1 = lines.pop(0)
    #print "<br>"
    #print "after fasta cleaning", "|".join(lines)
    #print "<br>"

    newSeq = []
    nCount = 0
    for l in lines:
        if len(l)==0:
            continue
        for c in l:
            if c not in "actgACTGNn":
                nCount +=1
            else:
                newSeq.append(c)
    seq = "".join(newSeq)

    msgs = []
    if len(seq)>MAXSEQLEN and db!="noGenome":
        msgs.append("<strong>Sorry, this tool cannot handle sequences longer than 1kbp</strong><br>Below you find the results for the first %d bp of your input sequence.<br>" % MAXSEQLEN)
        seq = seq[:MAXSEQLEN]
    if len(seq)>MAXSEQLEN_NOGENOME and db=="noGenome":
        msgs.append("<strong>Sorry, this tool cannot handle sequences longer than %d bp when specifying 'No Genome'.</strong><br>Below you find the results for the first %d bp of your input sequence.<br>" % (MAXSEQLEN_NOGENOME, MAXSEQLEN_NOGENOME))
        seq = seq[:MAXSEQLEN_NOGENOME]

    if nCount!=0:
        msgs.append("Sequence contained %d non-ACTGN letters. They were removed." % nCount)

    return seq, "<br>".join(msgs)

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M', 
        "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n"}
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def findPams (seq, pam, strand, startDict, endSet):
    """ return two values: dict with pos -> strand of PAM and set of end positions of PAMs
    Makes sure to return only values with at least GUIDELEN bp left (if strand "+") or to the
    right of the match (if strand "-")
    If the PAM is TTTN, then this is inversed: pos-strand matches must have at least GUIDELEN
    basepairs to the right, neg-strand matches must have at least GUIDELEN bp on their left
    >>> findPams("GGGGGGGGGGGGGGGGGGGGGGG", "NGG", "+", {}, set(), False)
    ({20: '+'}, set([23]))
    >>> findPams("CCAGCCCCCCCCCCCCCCCCCCC", "CCA", "-", {}, set(), False)
    ({0: '-'}, set([3]))
    >>> findPams("TTTNCCCCCCCCCCCCCCCCCTTTN", "TTTN", "+", {}, set(), True)
    ({0: '+'}, set([3]))
    >>> findPams("CCCCCCCCCCCCCCCCCCCCCAAAA", "NAA", "-", {}, set(), True
    ({}, set([]))
    >>> findPams("AAACCCCCCCCCCCCCCCCCCCCC", "NAA", "-", {}, set(), True)
    ({}, set([]))
    >>> findPams("CCCCCCCCCCCCCCCCCCCCCCCCCAA", "NAA", "-", {}, set(), True)
    ({}, set([]))

    """
    GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(pam)
    assert(cpf1Mode is not None)

    if cpf1Mode:
        maxPosPlus  = len(seq)-(GUIDELEN+len(pam))
        minPosMinus = GUIDELEN
    else:
        # -------------------
        #          OKOKOKOKOK
        minPosPlus  = GUIDELEN
        # -------------------
        # OKOKOKOKOK
        maxPosMinus = len(seq)-(GUIDELEN+len(pam))

    #print "new search", seq, pam, "<br>"
    for start in findPat(seq, pam):
        if cpf1Mode:
            # need enough flanking seq on one side
            #print "found", start,"<br>"
            if strand == "+" and start > maxPosPlus:
                continue
            if strand == "-" and start < minPosMinus:
                continue
        else:
            if strand=="+" and start < minPosPlus:
                continue
            if strand=="-" and start > maxPosMinus:
                continue

        #print "match", strand, start, end, "<br>"
        startDict[start] = strand
        end = start+len(pam)
        endSet.add(end)
    return startDict, endSet

def matchNuc(pat, nuc):
    " returns true if pat (single char) matches nuc (single char) "
    if pat in ["A", "C", "T", "G"] and pat==nuc:
        return True
    elif pat=="S" and nuc in ["G", "C"]:
        return True
    elif pat=="M" and nuc in ["A", "C"]:
        return True
    elif pat=="D" and nuc in ["AGT"]:
        return True
    elif pat=="W" and nuc in ["A", "T"]:
        return True
    elif pat=="K" and nuc in ["T", "G"]:
        return True
    elif pat=="R" and nuc in ["A", "G"]:
        return True
    elif pat=="Y" and nuc in ["C", "T"]:
        return True
    else:
        return False