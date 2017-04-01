
from constants import addGenePlasmidsAureus


def errAbort(msg):
    " print err msg and exit "
    raise Exception(msg)
    sys.exit(0)  # cgi must not exit with 1


def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M', 
        "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n"}
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def makeTempFile(prefix, suffix):
    " return a temporary file that is deleted upon exit, unless DEBUG is set "
    if DEBUG:
        fname = join("/tmp", prefix+suffix)
        fh = open(fname, "w")
    else:
        fh = tempfile.NamedTemporaryFile(prefix="primer3In", suffix=".txt")
    return fh

def setupPamInfo(pam):
    " modify a few globals based on the current pam "
    if pam=="TTTN":
        GUIDELEN = 23
        cpf1Mode = True
    elif pam=="NNGRRT":
        addGenePlasmids = addGenePlasmidsAureus
        GUIDELEN = 21
        cpf1Mode = False
    else:
        GUIDELEN = 20
        cpf1Mode = False
        addGenePlasmids = 'None'
    return GUIDELEN,cpf1Mode,addGenePlasmids
