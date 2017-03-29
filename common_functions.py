import hashlib, base64,string

from constants import addGenePlasmidsAureus

transTab = string.maketrans("-=/+_", "abcde")
def errAbort(msg):
    " print err msg and exit "
    raise Exception(msg)
    sys.exit(0)  # cgi must not exit with 1


def makeTempBase(seq, org, pam, batchName):
    "create the base name of temp files using a hash function and some prettyfication "
    hasher = hashlib.sha1(seq+org+pam+batchName)
    batchId = base64.urlsafe_b64encode(hasher.digest()[0:20]).translate(transTab)[:20]
    return batchId

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
