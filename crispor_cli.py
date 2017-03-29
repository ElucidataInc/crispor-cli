#!/usr/bin/env python2.7
# the tefor crispr tool
# can be run as a CGI or from the command line

# OOF scores are WRONG for Cpf1! -> where is the cut site?

# python std library
import subprocess, tempfile, optparse, logging, atexit, glob, shutil
import Cookie, time, sys, cgi, re, random, platform, os
import hashlib, base64, string, logging, operator, urllib, sqlite3, time
import traceback, json, pwd, pickle

from datetime import datetime
from collections import defaultdict, namedtuple
from os.path import join, isfile, basename, dirname, isdir, abspath
from StringIO import StringIO
from itertools import product

# try to load external dependencies
# we're going into great lengths to create a readable error message
needModules = set(["tabix", "twobitreader", "pandas", "matplotlib", "scipy"])
try:
    import tabix # if not found, install with 'pip install pytabix'
    needModules.remove("tabix") 
except:
    pass

try:
    import twobitreader # if not found, install with 'pip install twobitreader'
    needModules.remove("twobitreader")
except:
    pass
    
try:
    import pandas # required by doench2016 score. install with 'pip install pandas'
    needModules.remove("pandas")
    import scipy # required by doench2016 score. install with 'pip install pandas'
    needModules.remove("scipy")
    import matplotlib # required by doench2016 score. install with 'pip install matplotlib'
    needModules.remove("matplotlib")
    import numpy # required by doench2016 score. install with 'pip install numpy'
    needModules.remove("numpy")
except:
    pass

if len(needModules)!=0:
    print("Content-type: text/html\n")
    print("Python interpreter path: %s<p>" % sys.executable)
    print("These python modules were not found: %s<p>" % ",".join(needModules))
    print("To install all requirements in one line, run: pip install pytabix pandas twobitreader scipy matplotlib numpy<p>")
    sys.exit(0)

# our own eff scoring library
import crisporEffScores

# don't report print as an error
# pylint: disable=E1601

# optional module for Excel export as native .xls files
# install with 'apt-get install python-xlwt' or 'pip install xlwt'
xlwtLoaded = True
try:
    import xlwt
except:
    sys.stderr.write("crispor.py - warning - the python xlwt module is not available\n")
    xlwtLoaded = False

# optional module for mysql support
try:
    import MySQLdb
    mysqldbLoaded = True
except:
    mysqldbLoaded = False

# version of crispor
versionStr = "4.0"

# contact email
contactEmail='crispor@tefor.net'

# write debug output to stdout
DEBUG = False
#DEBUG = True

# use bowtie for off-target search?
useBowtie = False

# calculate the efficienc scores?
doEffScoring = True

# system-wide temporary directory
TEMPDIR = os.environ.get("TMPDIR", "/tmp")
#TEMPDIR = "/dev/shm/"

# prefix in html statements before the directories "image/", "style/" and "js/" 
HTMLPREFIX =  ""
# alternative directory on local disk where image/, style/ and js/ are located
HTMLDIR = "/usr/local/apache/htdocs/crispor/"

# directory of crispor.py
baseDir = dirname(__file__)

# filename of this script, usually crispor.py
myName = basename(__file__)

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

# directory for processed batches of offtargets ("cache" of bwa results)
batchDir = join(baseDir,"temp")

# use mysql or sqlite for the jobs?
jobQueueUseMysql = False

# the file where the sqlite job queue is stored
JOBQUEUEDB = join("/tmp/crisporJobs.db")

# alternatively: connection info for mysql
jobQueueMysqlConn = {"socket":None, "host":None, "user": None, "password" : None}

# directory for platform-independent scripts (e.g. Heng Li's perl SAM parser)
scriptDir = join(baseDir, "bin")

# directory for helper binaries (e.g. BWA)
binDir = abspath(join(baseDir, "bin", platform.system()))

# directory for genomes
genomesDir = join(baseDir, "genomes")

DEFAULTORG = 'hg19'
DEFAULTSEQ = 'cttcctttgtccccaatctgggcgcgcgccggcgccccctggcggcctaaggactcggcgcgccggaagtggccagggcgggggcgacctcggctcacagcgcgcccggctattctcgcagctcaccatgGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC'

# used if hg19 is not available
ALTORG = 'sacCer3'
ALTSEQ = 'ATTCTACTTTTCAACAATAATACATAAACatattggcttgtggtagCAACACTATCATGGTATCACTAACGTAAAAGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTTCAGAATATTTCGTACTTACACAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCT'

pamDesc = [ ('NGG','20bp-NGG - Cas9 S. Pyogenes'),
         ('TTTN','TTTN-23bp - Cpf1 F. novicida'),
         ('NGA','20bp-NGA - Cas9 S. Pyogenes mutant VQR'),
         ('NGCG','20bp-NGCG - Cas9 S. Pyogenes mutant VRER'),
         ('NNAGAA','20bp-NNAGAA - Cas9 S. Thermophilus'),
         ('NGGNG','20bp-NGGNG - Cas9 S. Thermophilus'),
         ('NNGRRT','21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus'),
         ('NNNNGMTT','20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis'),
         ('NNNNACA','20bp-NNNNACA - Cas9 Campylobacter jejuni'),
       ]

DEFAULTPAM = 'NGG'


# maximum size of an input sequence 
MAXSEQLEN = 1000
# maximum input size when specifying "no genome"
MAXSEQLEN_NOGENOME = 25000

# BWA: allow up to X mismatches
maxMMs=4

# maximum number of occurences in the genome to get flagged as repeats. 
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000/MAXOCC

# the length of the guide sequence, set by setupPamInfo
GUIDELEN=None

# the name of the currently processed batch, assigned only once 
# in readBatchParams and only for json-type batches
batchName = ""

# are we doing a Cpf1 run?
# this variable changes almost all processing and
# has to be set on program start, as soon as we know 
# the PAM we're running on
cpf1Mode=None


# Highly-sensitive mode (not for CLI mode): 
# MAXOCC is increased in processSubmission() and in the html UI if only one
# guide seq is run
# Also, the number of allowed mismatches is increased to 5 instead of 4
#HIGH_MAXOCC=600000
#HIGH_maxMMs=5

# minimum off-target score of standard off-targets (those that end with NGG)
# This should probably be based on the CFD score these days
# But for now, I'll let the user do the filtering
MINSCORE = 0.0

# minimum off-target score for alternative PAM off-targets
# There is not a lot of data to support this cutoff, but it seems
# reasonable to have at least some cutoff, as otherwise we would show
# NAG and NGA like NGG and the data shows clearly that the alternative
# PAMs are not recognized as well as the main NGG PAM.
# so for now, I just filter out very degenerative ones. the best solution
# would be to have a special penalty on the CFD score, but CFS does not 
# support non-NGG PAMs (is this actually true?)
ALTPAMMINSCORE = 1.0

# for some PAMs, we allow other alternative motifs when searching for offtargets
# MIT and eCrisp do that, they use the motif NGG + NAG, we add one more, based on the
# on the guideSeq results in Tsai et al, Nat Biot 2014
# The NGA -> NGG rule was described by Kleinstiver...Young 2015 "Improved Cas9 Specificity..."
# NNGTRRT rule for S. aureus is in the new protocol "SaCas9 User manual"
# ! the length of the alternate PAM has to be the same as the original PAM!
offtargetPams = {
"NGG" : ["NAG","NGA"],
"NGA" : ["NGG"],
"NNGRRT" : ["NNGRRN"]
}

# global flag to indicate if we're run from command line or as a CGI
commandLineMode = False

# names/order of efficiency scores to show in UI
scoreNames = ["fusi", "crisprScan"]
allScoreNames = ["fusi", "chariRank", "ssc", "doench", "wang", "crisprScan", "oof", "housden", "proxGc"]

# how many digits shall we show for each score? default is 0
scoreDigits = {
    "ssc" : 1,
}

# List of AddGene plasmids, their long and short names:
addGenePlasmids = [
("43860", ("MLM3636 (Joung lab)", "MLM3636")),
("49330", ("pAc-sgRNA-Cas9 (Liu lab)", "pAcsgRnaCas9")),
("42230", ("pX330-U6-Chimeric_BB-CBh-hSpCas9 (Zhang lab) + derivatives", "pX330")),
("52961", ("lentiCRISPR v2 (Zhang lab)", "lentiCrispr")),
]

addGenePlasmidsAureus = [
("61591", ("pX601-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA;U6::BsaI-sgRNA (Zhang lab)", "pX601")),
("61592", ("pX600-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA (Zhang lab)", "pX600")),
("61593", ("pX602-AAV-TBG::NLS-SaCas9-NLS-HA-OLLAS-bGHpA;U6::BsaI-sgRNA (Zhang lab)", "pX602")),
("65779", ("VVT1 (Joung lab)", "VVT1"))
]

# list of AddGene primer 5' and 3' extensions, one for each AddGene plasmid
# format: prefixFw, prefixRw, restriction enzyme, link to protocol
addGenePlasmidInfo = {
"43860" : ("ACACC", "AAAAC", "BsmBI", "https://www.addgene.org/static/data/plasmids/43/43860/43860-attachment_T35tt6ebKxov.pdf"),
"49330" : ("TTC", "AAC", "Bsp QI", "http://bio.biologists.org/content/3/1/42#sec-9"),
"42230" : ("CACC", "AAAC", "Bbs1", "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf"),
"52961" : ("CACC", "AAAC", "Bbs1", "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf"),
"61591" : ("CACC", "AAAC", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61591/61591-attachment_it03kn5x5O6E.pdf"),
"61592" : ("CACC", "AAAC", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf"),
"61593" : ("CACC", "AAAC", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf"),
"65779": ("CACC", "AAAC", "BsmBI", "https://www.addgene.org/static/data/plasmids/65/65779/65779-attachment_G8oNyvV6pA78.pdf")
}

# Restriction enzyme supplier codes
rebaseSuppliers = {
"B":"Life Technologies",
"C":"Minotech",
"E":"Agilent",
"I":"SibEnzyme",
"J":"Nippon Gene",
"K":"Takara",
"M":"Roche",
"N":"NEB",
"O":"Toyobo",
"Q":"Molecular Biology Resources",
"R":"Promega",
"S":"Sigma",
"V":"Vivantis",
"X":"EURx",
"Y":"SinaClon BioScience"
}

# labels and descriptions of eff. scores
scoreDescs = {
    "doench" : ("Doench '14", "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>"),
    "ssc" : ("Xu", "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>"),
    "crisprScan" : ["Moreno-Mateos", "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score."],
    "wang" : ("Wang", "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>"),
    "chariRank" : ("Chari", "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>"),
    "fusi" : ("Doench '16", "Previously called the 'Fusi' score. Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score."),
    "housden" : ("Housden", "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>"),
    "oof" : ("Out-of-Frame", "Range: 0-100. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al.</a>. Click the score to show the most likely deletions for this guide.")
}

# the headers for the guide and offtarget output files
guideHeaders = ["guideId", "guideSeq", "mitSpecScore", "cfdSpecScore", "offtargetCount", "guideGenomeMatchGeneLocus"]
offtargetHeaders = ["guideId", "guideSeq", "offtargetSeq", "mismatchPos", "mismatchCount", "mitOfftargetScore", "cfdOfftargetScore", "chrom", "start", "end", "strand", "locusDesc"]

# a file crispor.conf in the directory of the script allows to override any global variable
myDir = dirname(__file__)
confPath =join(myDir, "crispor.conf")
if isfile(confPath):
    exec(open(confPath))

cgiParams = None

# ====== END GLOBALS ============

def setupPamInfo(pam):
    " modify a few globals based on the current pam "
    global GUIDELEN
    global cpf1Mode
    global addGenePlasmids
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


# ==== CLASSES =====
class JobQueue:
    """
    simple job queue, using a db table as a backend
    jobs have different types and status. status can be updated while they run
    job running times are kept and old job info is kept in a separate table
    
    >>> q = JobQueue()
    >>> q.clearJobs()
    >>> q.waitCount()
    0
    >>> q.addJob("search", "abc123", "myParams")
    True

    only one job per jobId
    >>> q.addJob("search", "abc123", "myParams")
    False
    >>> q.waitCount()
    1
    >>> q.getStatus("abc123")
    u'Waiting'
    >>> q.startStep("abc123", "bwa", "Alignment with BWA")
    >>> q.getStatus("abc123")
    u'Alignment with BWA'
    >>> jobType, jobId, paramStr = q.popJob()

    >>> q.waitCount()
    0
    >>> q.jobDone("abc123")
    >>> q.waitCount()
    0

    can't pop from an empty queue
    #>>> q.popJob()
    #(None, None, None)
    #>>> os.system("rm /tmp/tempCrisporTest.db")
    #0
    """

    _queueDef = (
    'CREATE TABLE IF NOT EXISTS %s '
    '('
    '  jobType text,' # either "index" or "search"
    '  jobId text %s,' # unique identifier
    '  paramStr text,' # parameters for jobs, like db, options, etc. 
    '  isRunning int DEFAULT 0,' # indicates steps have started, done jobs are moved to doneJobs table
    '  stepName text,' # currently step, internal step name for timings
    '  stepLabel text,' # current step, human-readable status of job, for UI
    '  lastUpdate float,' # time of last update
    '  stepTimes text,' # comma-sep list of whole msecs, one per step
    '  startTime text' # date+time when job was put into queue
    ')')

    def __init__(self):
        " no inheritance needed here "
        if jobQueueUseMysql:
            self.openMysql()
        else:
            self.openSqlite(JOBQUEUEDB)

    def openSqlite(self, dbName):
        self.dbName = dbName
        self.conn = sqlite3.Connection(dbName)
        try:
            self.conn.execute(self._queueDef % ("queue", "PRIMARY KEY"))
        except sqlite3.OperationalError:
            errAbort("cannot open the file %s" % JOBQUEUEDB)
        self.conn.execute(self._queueDef % ("doneJobs", ""))
        self.conn.commit()
        self._chmodJobDb()

    def _chmodJobDb(self):
        # umask is not respected by sqlite, bug http://www.mail-archive.com/sqlite-users@sqlite.org/msg59080.html
        if not jobQueueUseMysql:
            try:
                os.chmod(JOBQUEUEDB, 0o666)
            except OSError:
                pass

    def openMysql(self):
        db = MySQLdb.connect(**jobQueueMysql)
        self.conn = db.cursor()

    def addJob(self, jobType, jobId, paramStr):
        " create a new job, returns False if not successful  "
        self._chmodJobDb()
        sql = 'INSERT INTO queue (jobType, jobId, isRunning, lastUpdate, ' \
            'stepTimes, paramStr, stepName, stepLabel, startTime) VALUES (?, ?, ?, ?, ?, ?, ?, ?, datetime("now"))'
        now = "%.3f" % time.time()
        try:
            self.conn.execute(sql, (jobType, jobId, 0, now, "", paramStr, "wait", "Waiting"))
            self.conn.commit()
            return True
        except sqlite3.IntegrityError:
            return False
        except sqlite3.OperationalError:
            errAbort("Cannot open DB file %s. Please contact %s" % (self.dbName, contactEmail))

    def getStatus(self, jobId):
        " return current job status label or None if job is not in queue"
        sql = 'SELECT stepLabel FROM queue WHERE jobId=?'
        try:
            status = self.conn.execute(sql, (jobId,)).next()[0]
        except StopIteration:
            status = None
        return status

    def dump(self):
        " for debugging, write the whole queue table to stdout "
        sql = 'SELECT * FROM queue'
        for row in self.conn.execute(sql):
            print "\t".join([str(x) for x in row])

    def jobInfo(self, jobId, isDone=False):
        " for debugging, return all job info as a tuple "
        if isDone:
            sql = 'SELECT * FROM doneJobs WHERE jobId=?'
        else:
            sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = self.conn.execute(sql, (jobId,)).next()
        except StopIteration:
            return []
        return row

    def startStep(self, jobId, newName, newLabel):
        " start a new step. Update lastUpdate, status and stepTime "
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT lastUpdate, stepTimes, stepName FROM queue WHERE jobId=?'
        lastTime, timeStr, lastStep = self.conn.execute(sql, (jobId,)).next()
        lastTime = float(lastTime)

        # append a string in format "stepName:milliSecs" to the timeStr
        now = time.time()
        timeDiff = "%d" % int((1000.0*(now - lastTime)))
        newTimeStr = timeStr+"%s=%s" % (lastStep, timeDiff)+","

        sql = 'UPDATE queue SET lastUpdate=?, stepName=?, stepLabel=?, stepTimes=?, isRunning=? WHERE jobId=?'
        self.conn.execute(sql, (now, newName, newLabel, newTimeStr, 1, jobId))

        tryCount = 0
        while tryCount < 10:
            try:
                self.conn.commit()
                break
            except sqlite3.OperationalError:
                time.sleep(3)
                tryCount += 1

        if tryCount >= 10:
            raise Exception("Database locked for a long time")

    def jobDone(self, jobId):
        " remove the job from the queue and add it to the queue log"
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = self.conn.execute(sql, (jobId,)).next()
        except StopIteration:
            # return if the job has already been removed
            logging.warn("jobDone - jobs %s has been removed already" % jobId)
            return

        sql = 'DELETE FROM queue WHERE jobId=?'
        self.conn.execute(sql, (jobId,))

        sql = 'INSERT INTO doneJobs VALUES (?,?,?,?,?,?,?,?,?)'
        self.conn.execute(sql, row)
        self.conn.commit()

    def waitCount(self):
        " return number of waiting jobs "
        sql = 'SELECT count(*) FROM queue WHERE isRunning=0'
        count = self.conn.execute(sql).next()[0]
        return count

    def popJob(self):
        " return (jobType, jobId, params) of first waiting job and set it to running state "
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT jobType, jobId, paramStr FROM queue WHERE isRunning=0 ORDER BY lastUpdate LIMIT 1'
        try:
            jobType, jobId, paramStr = self.conn.execute(sql).next()
        except StopIteration:
            self.conn.commit() # unlock db
            return None, None, None

        sql = 'UPDATE queue SET isRunning=1 where jobId=?'
        self.conn.execute(sql, (jobId,))
        self.conn.commit() # unlock db
        return jobType, jobId, paramStr

    def clearJobs(self):
        " clear the job table, removing running jobs, too "
        self.conn.execute("DELETE from queue")
        #self.conn.execute("DROP TABLE queue")
        self.conn.commit()

    def close(self):
        " "
        self.conn.close()

# ====== FUNCTIONS =====
contentLineDone = False

def errAbort(msg):
    " print err msg and exit "
    if commandLineMode:
        raise Exception(msg)

    if not contentLineDone:
        print "Content-type: text/html\n"

    print('<div style="float:left; text-align:left; width: 800px">')
    print("<strong>Error:</strong> ")
    print(msg+"<p>")
    print("If you think this is a bug or you have any other suggestions, please do not hesitate to email %s" % contactEmail)
    print('</div>')
    sys.exit(0)  # cgi must not exit with 1

# allow only dashes, digits, characters, underscores and colons in the CGI parameters
# and +
notOkChars = re.compile(r'[^a-zA-Z0-9-_:+.]')

def checkVal(key, str):
    """ remove special characters from input string, to protect against injection attacks """
    if len(str) > 10000:
	errAbort("input parameter %s is too long" % key)
    if notOkChars.search(str):
	errAbort("input parameter %s contains an invalid character" % key)
    return str

def cgiGetParams():
    " get CGI parameters and return as dict "
    form = cgi.FieldStorage()
    global cgiParams
    cgiParams = {}

    # parameters are:
    #"pamId", "batchId", "pam", "seq", "org", "download", "sortBy", "format", "ajax
    for key in form.keys():
        val = form.getfirst(key)
	if val!=None:
            # "seq" is cleaned by cleanSeq later
            val = urllib.unquote(val)
            if key not in ["seq", "name"]:
                checkVal(key, val)
            cgiParams[key] = val

    if "pam" in cgiParams:
        if len(set(cgiParams["pam"])-set("ACTGNMKRY"))!=0:
            errAbort("Illegal character in PAM-sequence. Only ACTGMKRY and N allowed.")
    return cgiParams

transTab = string.maketrans("-=/+_", "abcde")

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

def saveSeqOrgPamToCookies(seq, org, pam):
    " create a cookie with seq, org and pam and print it"
    cookies=Cookie.SimpleCookie()
    expires = 365 * 24 * 60 * 60
    cookies['lastseq'] = seq
    cookies['lastseq']['expires'] = expires
    cookies['lastorg'] = org
    cookies['lastorg']['expires'] = expires
    cookies['lastpam'] = pam
    cookies['lastpam']['expires'] = expires
    print cookies

def debug(msg):
    if commandLineMode:
        logging.debug(msg)
    elif DEBUG:
        print msg
        print "<br>"

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

def gcContent(seq):
    " return GC content as a float "
    c = 0
    for x in seq:
        if x in ["G","C"]:
            c+= 1
    return (float(c)/len(seq))
            
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

def rulerString(maxLen):
    " return line with positions every 10 chars "
    texts = []
    for i in range(0, maxLen, 10):
        numStr = str(i)
        texts.append(numStr)
        spacer = "".join([" "]*(10-len(numStr)))
        texts.append(spacer)
    return "".join(texts)

def varDictToHtml(varDict, seq):
    " make a list of one html string per position in the sequence "
    if varDict is None:
        return None

    varHtmls = []
    for i in range(0, len(seq)):
        if not i in varDict:
            varHtmls.append(".")
        else:
            varHooverLines = []
            showStar = False # show a star if change is non-simple SNP
            varInfos = varDict[i]
            for chrom, pos, refAll, altAll, infoDict in varInfos:
                varHooverLines.append("%s &rarr; %s<br>" % (refAll, altAll))
                if "freq" in infoDict:
                    varHooverLines.append("&nbsp;<b>Freq:</b> %s<br>" % infoDict["freq"])
                #if "dbg" in infoDict:
                    #varHooverLines.append("%s<br>" % infoDict["dbg"])
                if "varId" in infoDict:
                    varHooverLines.append("&nbsp;<b>ID:</b> %s<br>" % infoDict["varId"])
                if "ExAC" in varDict["label"]:
                        endPos = int(pos)+len(refAll)
                        varHooverLines.append('&nbsp;<a target=_blank href="http://exac.broadinstitute.org/region/%s-%s-%d">ExAC Browser</a><br>' % (chrom, pos, endPos))

                if len(refAll)!=1 or len(altAll)!=1:
                    showStar = True

            if len(varInfos)!=1:
                showStar = True

            varDesc = "".join(varHooverLines)

            if showStar:
                dispChar = "*"
            else:
                dispChar = altAll
            varHtmls.append("<u class='tooltipsterInteract' title='%s'>%s</u>" % (varDesc, dispChar))
    return varHtmls

def showSeqAndPams(seq, startDict, pam, guideScores, varHtmls, varDbs, varDb, minFreq):
    " show the sequence and the PAM sites underneath in a sequence viewer "
    lines, maxY = distrOnLines(seq.upper(), startDict, len(pam))

    print "<div class='substep'>"
    print '<a id="seqStart"></a>'
    print "Found %d possible guide sequences in input (%d bp). Click on a PAM %s match to show its guide sequence.<br>" % (len(guideScores), len(seq), pam)

    if not cpf1Mode:
        print "Shown below are the PAM site and the expected cleavage position located -3bp 5' of the PAM site.<br>"
        print '''Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.<br>'''
    else:
        print ""

    if varDb is not None:
        print("""<form style="display:inline" id="paramForm" action="%s" method="GET">""" % basename(__file__))
        print ("Variant database:")
        varDbList = [(b,c) for a,b,c,d in varDbs] # only keep fname+label
        printDropDown("varDb", varDbList, varDb)

        if minFreq==0.0:
            minFreq="0.0"
        else:
            minFreq = str(minFreq)

        # pull out the hasAF field for this varDb
        varDbHasAF = False
        for shortLabel, fname, desc, hasAF in varDbs:
            if fname==varDb:
                varDbHasAF = hasAF
                break

        if varDbHasAF:
            print(""" Min. frequency: """)
            print("""<input type="text" name="minFreq" size="8" value="%s">""" % minFreq)
        print("""<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">""")


    print "</div>"
    print '''<div style="text-align: left; overflow-x:scroll; width:100%; background:#DDDDDD; border-style: solid; border-width: 1px">'''

    print '<pre style="font-size: 80%; display:inline; line-height: 0.95em; text-align:left">'+rulerString(len(seq))

    if varHtmls is not None:
        print "".join(varHtmls)
    print seq

    for y in range(0, maxY+1):
        texts = []
        lastEnd = 0
        for start, end, name, strand, pamId  in lines[y]:
            spacer = "".join([" "]*((start-lastEnd)))
            lastEnd = end
            texts.append(spacer)
            score = guideScores[pamId]
            color = scoreToColor(score)

            texts.append('''<a style="text-shadow: 1px 1px 1px #bbb; color: %s" id="list%s" href="#%s">''' % (color, pamId,pamId))
            texts.append(name)
            texts.append("</a>")
        print(u''.join(texts).encode("utf8"))
    print("</pre><br>")

    print '''</div>'''

    if cpf1Mode:
        print('<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">Cpf1 has a staggered site: cleavage occurs after the 18th base on the non-targeted strand which has the TTTN PAM motif (indicate by "\\" in the schema above). Cleavage occurs after the 23rd base on the targeted strand which has the AAAN motif (indicated by "/" in the schema above). See <a target=_blank href="http://www.sciencedirect.com/science/article/pii/S0092867415012003">Zetsche et al 2015</a>, in particular <a target=_blank href="http://www.sciencedirect.com/science?_ob=MiamiCaptionURL&_method=retrieve&_eid=1-s2.0-S0092867415012003&_image=1-s2.0-S0092867415012003-gr3.jpg&_cid=272196&_explode=defaultEXP_LIST&_idxType=defaultREF_WORK_INDEX_TYPE&_alpha=defaultALPHA&_ba=&_rdoc=1&_fmt=FULL&_issn=00928674&_pii=S0092867415012003&md5=11771263f3e390e444320cacbcfae323">Fig 3</a>.</div>')
    
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

def flankSeqIter(seq, startDict, pamLen, doFilterNs):
    """ given a seq and dictionary of pos -> strand and the length of the pamSite
    yield tuples of (name, pamStart, guideStart, strand, flankSeq, pamSeq)

    if doFilterNs is set, will not return any sequences that contain an N character
    """
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

def makeBrowserLink(dbInfo, pos, text, title, cssClasses):
    " return link to genome browser (ucsc or ensembl) at pos, with given text "
    if dbInfo.server.startswith("Ensembl"):
        baseUrl = "www.ensembl.org"
        urlLabel = "Ensembl"
        if dbInfo.server=="EnsemblPlants":
            baseUrl = "plants.ensembl.org"
        elif dbInfo.server=="EnsemblMetazoa":
            baseUrl = "metazoa.ensembl.org"
        elif dbInfo.server=="EnsemblProtists":
            baseUrl = "protists.ensembl.org"
        org = dbInfo.scientificName.replace(" ", "_")
        url = "http://%s/%s/Location/View?r=%s" % (baseUrl, org, pos)
    elif dbInfo.server=="ucsc":
        urlLabel = "UCSC"
        if pos[0].isdigit():
            pos = "chr"+pos
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (dbInfo.name, pos)
    # some limited support for gbrowse
    elif dbInfo.server.startswith("http://"):
        urlLabel = "GBrowse"
        chrom, start, end, strand = parsePos(pos)
        start = start+1
        url = "%s/?name=%s:%d..%d" % (dbInfo.server, chrom, start, end)
    else:
        #return "unknown genome browser server %s, please email services@tefor.net" % dbInfo.server
        urlLabel = None
        url = "javascript:void(0)"

    classStr = ""
    if len(cssClasses)!=0:
        classStr = ' class="%s"' % (" ".join(cssClasses))

    if title is None:
        if urlLabel != None:
            title = "Link to %s Genome Browser" % urlLabel
        else:
            title = "No Genome Browser link available yet for this organism"
    return '''<a title="%s"%s target="_blank" href="%s">%s</a>''' % (title, classStr, url, text)

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

def parseChromSizes(genome):
    " return chrom sizes as dict chrom -> size "
    genomeDir = genomesDir # make local
    sizeFname = "%(genomeDir)s/%(genome)s/%(genome)s.sizes" % locals()
    ret = {}
    for line in open(sizeFname).read().splitlines():
        fields = line.split()
        chrom, size = fields[:2]
        ret[chrom] = int(size)
    return ret

def extendAndGetSeq(db, chrom, start, end, strand, flank=100):
    """ extend (start, end) by flank and get sequence for it using twoBitTwoFa.
    Return None if not possible to extend.
    #>>> extendAndGetSeq("hg19", "chr21", 10000000, 10000005, "+", flank=3)
    #'AAGGAATGTAG'
    """
    chromSizes = parseChromSizes(db)
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

def pamStartToGuideRange(startPos, strand, pamLen):
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

def htmlHelp(text):
    " show help text with tooltip or modal dialog "
    className = "tooltipster"
    if "href" in text:
        className = "tooltipsterInteract"

    print '''<img style="height:1.1em; width:1.0em" src="%simage/info-small.png" class="help %s" title="%s" />''' % (HTMLPREFIX, className, text)

def htmlWarn(text):
    " show help text with tooltip "
    print '''<img style="height:1.1em; width:1.0em" src="%simage/warning-32.png" class="help tooltipster" title="%s" />''' % (HTMLPREFIX, text)

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

    pamSeqs = list(flankSeqIter(seq, startDict, len(pamPat), True))

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

def printDownloadTableLinks(batchId):
    print '<div style="text-align:right">'
    print '<small>'
    print "Download tables: "
    print '<a href="crispor.py?batchId=%s&download=guides&format=xls">Guides</a>&nbsp;' % batchId
    print '<a href="crispor.py?batchId=%s&download=offtargets&format=xls">Off-targets</a>' % batchId
    print '</small>'
    print '</div>'

def hasGeneModels(org):
    " return true if this organism has gene model information "
    geneFname = join(genomesDir, org, org+".segments.bed")
    return isfile(geneFname)

def printTableHead(batchId, chrom, org):
    " print guide score table description and columns "
    # one row per guide sequence
    if not cpf1Mode:
        print '''<div class='substep'>Ranked by default from highest to lowest specificity score (<a target='_blank' href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al., Nat Biot 2013</a>) as on <a href="http://crispr.mit.org">http://crispr.mit.org</a>. Click on a column title to rank by a different score.<br>'''
        print('''
        <b>Our recommendation:</b> Use Fusi for in-vivo (U6) transcribed guides, Moreno-Mateos for in-vitro (T7) guides injected into Zebrafish/Mouse oocytes.<br> See our <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2">CRISPOR paper in Gen Biol 2016</a>, figures 4 and 5.<br>''')
        #print('''References for scores:
        #<a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench/Fusi 2016</a>,
        #<a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos</a>''')

        #print(''',
        #<a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari</a>,
        #<a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu</a>,
        #<a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench</a>,
        #<a target='_blank' href='http://http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang</a>,
        #<a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden</a>,
        #<a target='_blank' href='http://www.cell.com/cell-reports/abstract/S2211-1247%2814%2900827-4'>Prox. GC</a>,
        #<a target='_blank' href='https://mcb.berkeley.edu/labs/meyer/publicationpdfs/959.full.pdf'>-GG</a>,
        #<a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Out-of-Frame</a>''')
        print('</div>')

    printDownloadTableLinks(batchId)

    print """
    <script type="text/javascript">
    function onlyExons() {
        if ($("#onlyExonBox").prop("checked")) { 
            $(".otMore").show();
            $(".otMoreLink").hide();
            $(".otLessLink").hide();
            $(".notExon").hide();
            }
        else {
            if ($("#onlySameChromBox").prop("checked")) {
                $(".notExon:not(.diffChrom)").show();
            }
            else {
                $(".notExon").show();
                $(".otMoreLink").show();
                $(".otMore").hide();
            }
        }
    }
    function onlySameChrom() {
        if ($("#onlySameChromBox").prop("checked"))
            { 
            $(".otMore").show();
            $(".otMoreLink").hide();
            $(".otLessLink").hide();
            $(".diffChrom").hide();
            }
        else {
            if ($("#onlyExonBox").prop("checked")) {
                $(".diffChrom:not(.notExon)").show();
            }
            else {
                $(".diffChrom").show();
                $(".otMoreLink").show();
                $(".otMore").hide();
            }
        }
    }

    function showAllOts(classId) {
        $("#"+classId).show();
        $("#"+classId+"MoreLink").hide();
        $("#"+classId+"LessLink").show();
    }
    function showLessOts(classId) {
        $("#"+classId).hide();
        $("#"+classId+"MoreLink").show();
        $("#"+classId+"LessLink").hide();
    }
    </script>
    """

    if not cpf1Mode:
        print '<table id="otTable" style="background:white;table-layout:fixed; overflow:scroll; width:100%">'
        #if len(scoreNames)!=2:
            #print '<colgroup span="1"><col></colgroup>'
            #print '<colgroup span="1"><col></colgroup>'
            #print '<colgroup span="1"><col></colgroup>'
            #print '<colgroup span="8"><col><col><col><col><col><col><col><col><col></colgroup>'
            #print '<colgroup span="1"><col></colgroup>'
            #print '<colgroup span="1"><col></colgroup>'
            #print '<colgroup span="1"><col></colgroup>'
    else:
        print '<table id="otTable" style="background:white;table-layout:fixed; overflow:scroll; width:100%">'

    print '<thead>'
    print '<tr style="border-bottom:none; border-left:5px solid black; background-color:#F0F0F0">'
    
    print '<th style="width:80px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=pos" class="tooltipster" title="Click to sort the table by the position of the PAM site">Position/<br>Strand</a>' % batchId
    htmlHelp("You can click on the links in this column to highlight the <br>PAM site in the sequence viewer at the top of the page.")
    print '</th>'

    print '<th style="width:180px; border-bottom:none">Guide Sequence + <i>PAM</i><br>Restriction Enzymes'

    if not cpf1Mode:
        print '<th style="width:70px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=spec" class="tooltipster" title="Click to sort the table by specificity score">Specificity Score</a>' % batchId
        htmlHelp("The higher the specificity score, the lower are off-target effects in the genome.<br>The specificity score ranges from 0-100 and measures the uniqueness of a guide in the genome. See <a href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al. Nat Biotech 2013</a>. We recommend values &gt;50, where possible.")
        print "</th>"

    if not cpf1Mode:
        width = "230px"

        if len(scoreNames)==2:
           print '<th style="width:150px; border-bottom:none" colspan="%d">Predicted Efficiency' % (len(scoreNames))
        else:
           print '<th style="width:230px; border-bottom:none" colspan="%d">Predicted Efficiency' % (len(scoreNames)-1) # -1 because proxGc is in scoreNames but has no col

        htmlHelp("The higher the efficiency score, the more likely is cleavage at this position. For details on the scores, mouseover their titles below.")

        if cgiParams.get("showAllScores", "0")=="0":
            print("""<a style="font-size:12px" href="%s" class="tooltipsterInteract" title="By default, only the two most relevant scores are shown, based on our study <a href='http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2'>Haeussler et al. 2016</a>. Click this link to show all efficiency scores.">Show all scores</a>""" % cgiGetSelfUrl({"showAllScores":"1"}, anchor="otTable"))
            scoreDescs["crisprScan"][0] = "Mor.-Mateos"
        else:
            print("""<a style="font-size:12px" href="%s" class="tooltipsterInteract" title="Show only the two main scores">Main scores</a>""" % cgiGetSelfUrl({"showAllScores":None}, anchor="otTable"))

    #print '<th style="width:50">Prox. GC'
    #htmlHelp("")
    print '</th>'

    if len(scoreNames)==2:
        oofWidth=100
        oofName="Out-of-Frame score"
        oofDesc = "Click on score to show micro-homology"
    else:
        oofWidth=45
        oofName="Out-of- Frame"
        oofDesc = "Click score for details"

    print '<th style="width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=oofScore" class="tooltipster" title="Click to sort the table by Out-of-Frame score">%s</a>' % (oofWidth, batchId, oofName)

    htmlHelp(scoreDescs["oof"][1])
    print "<br><br><small>%s</small>" % oofDesc
    print '</th>'

    print '<th style="width:117px; border-bottom:none">Off-targets for <br>0-1-2-3-4 mismatches<br><span style="color:grey">+ next to PAM </span>'
    htmlHelp("For each number of mismatches, the number of off-targets is indicated.<br>Example: 1-3-20-50-60 means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, <br>20 off-targets with 2 mismatches, etc.<br>Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.<br>Shown in grey are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.")
    #print "</th>"

    #print '<th style="width:120">Off-targets with no mismatches next to PAM</i>'
    print "</th>"
    print '<th style="width:*; border-bottom:none">Genome Browser links to matches sorted by CFD off-target score'
    htmlHelp("For each off-target the number of mismatches is indicated and linked to a genome browser. <br>Matches are ranked by CFD off-target score (see Doench 2016 et al) from most to least likely.<br>Matches can be filtered to show only off-targets in exons or on the same chromosome as the input sequence.")

    print '<br><small>'
    #print '<form id="filter-form" method="get" action="crispor.py#otTable">'
    print '<input type="hidden" name="batchId" value="%s">' % batchId

    if hasGeneModels(org):
        print '''<input type="checkbox" id="onlyExonBox" onchange="onlyExons()">exons only'''
    else:
        print '<small style="color:grey">No exons.</small>'

    if chrom!="":
        if chrom[0].isdigit():
            chrom = "chrom "+chrom
        print '''<input type="checkbox" id="onlySameChromBox" onchange="onlySameChrom()">%s only''' % chrom
    else:
        print '<small style="color:grey">&nbsp;No match, no chrom filter</small>'
    # a hidden submit button
    # print '<input  type="submit" name="submit" value="submit">'
    #print '<input  type="submit" name="submit" value="1" style="position: absolute; height: 0px; width: 0px; border: none; padding: 0px;" hidefocus="true" tabindex="-1"/>'
    #print '</form></small>'
    print "</small>"
    print "</th>"
    print "</tr>"

    # subheaders
    print '<tr style="border-top:none; border-left: solid black 5px; background-color:#F0F0F0">'

    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'

    if not cpf1Mode:
        print '<th style="border-top:none"></th>'

    if not cpf1Mode:
        for scoreName in scoreNames:
            if scoreName in ["oof", "proxGc"]:
                continue
            scoreLabel, scoreDesc = scoreDescs[scoreName]
            print '<th style="width: 10px; border: none; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>' % (scoreDesc, batchId, scoreName, scoreLabel)

    if not cpf1Mode and "proxGc" in scoreNames:
        # the ProxGC score comes next
        print '''<th style="border: none; border-top:none; border-right: none; border-left:none" class="rotate">'''
        print '''<div><span style="border-bottom:none">'''
        print '''<a title="This column shows two heuristics based on observations rather than computational models: <a href='http://www.cell.com/cell-reports/abstract/S2211-1247%2814%2900827-4'>Ren et al</a> 2014 obtained the highest cleavage in Drosophila when the final 6bp contained &gt;= 4 GCs, based on data from 39 guides. <a href='http://www.genetics.org/content/early/2015/02/18/genetics.115.175166.abstract'>Farboud et al.</a> obtained the highest cleavage in C. elegans for the 10 guides that ended with -GG, out of the 50 guides they tested.<br>The column contains + if the final GC count is &gt;= 4 and GG if the guide ends with GG." href="crispor.py?batchId=%s&sortBy=finalGc6" class="tooltipsterInteract">Prox GC</span></div></th>''' % (batchId)

    # these are empty cells to fill up the row and avoid white space
    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'
    print "</tr>"
    print '</thead>'

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

def makeOtBrowserLinks(otData, chrom, dbInfo, pamId):
    " return a list with the html texts of the offtarget links "
    links = []

    i = 0
    for otSeq, score, cfdScore, editDist, pos, gene, alnHtml in otData:
        cssClasses = ["tooltipster"]
        if not gene.startswith("exon:"):
            cssClasses.append("notExon")
        if pos.split(":")[0]!=chrom:
            cssClasses.append("diffChrom")

        classStr =  ""
        if len(cssClasses)!=0:
            classStr = ' class="%s"' % " ".join(cssClasses)

        link = makeBrowserLink(dbInfo, pos, gene, alnHtml, cssClasses)
        editDist = str(editDist)
        links.append( '''<div%(classStr)s>%(editDist)s:%(link)s</div>''' % locals() )

        #i+=1
        #if i>=3:
            #break

    #if not showAll and len(otData)>3:
         #print '''... <br>&nbsp;&nbsp;&nbsp;<a href="crispor.py?batchId=%s&showAll=1">- show all offtargets</a>''' % batchId
    return links

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

def printNoEffScoreFoundWarn(effScoresCount):
    if effScoresCount==0 and not cpf1Mode:
        print('<div style="text-align:left"><strong>Note:</strong> No guide could be scored for efficiency.')
        print("This happens when the input sequence is shorter than 100bp and there is no genome available to extend it. Please add flanking 50bp on both sides of the input sequence and submit this new, longer sequence.</div><br>")


def showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, chrom, varHtmls):
    " shows table of all PAM motif matches "
    print "<br><div class='title'>Predicted guide sequences for PAMs</div>" 

    global scoreNames
    if (cgiParams.get("showAllScores", "0")=="1"):
        scoreNames = allScoreNames

    showPamWarning(pam)
    showNoGenomeWarning(dbInfo)
    printTableHead(batchId, chrom, org)

    count = 0
    effScoresCount = 0
    showProxGcCol = ("proxGc" in scoreNames)
    
    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, pamStart, guideStart, strand, pamId, guideSeq, \
            pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow

        color = scoreToColor(guideScore)
        print '<tr id="%s" style="border-left: 5px solid %s">' % (pamId, color)

        # position and strand
        #print '<td id="%s">' % pamId
        print '<td>'
        print '<a href="#list%s">' % (pamId)
        print str(pamStart+1)+" /"
        if strand=="+":
            print 'fw'
        else:
            print 'rev'
        print '</a>'
        print "</td>"

        # sequence with variants and PCR primer link
        print "<td>"
        print "<small>"

        # guide sequence + PAM sequence
        if cpf1Mode:
            fullGuideHtml = "<tt><i>"+pamSeq+"</i> " + guideSeq+"</tt>"
            spacePos = len(pamSeq)
        else:
            fullGuideHtml = "<tt>"+guideSeq + " <i>" + pamSeq+"</i></tt>"
            spacePos = len(guideSeq)

        # variant-string
        if varHtmls is not None:
            varFound = False
            varStrs = []
            for i in range(guideStart, guideStart+len(guideSeq)+len(pamSeq)):
                html = varHtmls[i]
                if html!=".":
                    varFound = True
                if (i-guideStart)==spacePos:
                    varStrs.append("&nbsp;")
                varStrs.append(html)
            print("<tt>%s</tt><br>" % ("".join(varStrs)))

        print fullGuideHtml
        print "<br>"

        scriptName = basename(__file__)
        if otData!=None and subOptMatchCount <= MAXOCC:
            print('<a href="%s?batchId=%s&pamId=%s&pam=%s" target="_blank">PCR primers</a>' % (scriptName, batchId, urllib.quote(str(pamId)), pam) )

        if gcContent(guideSeq)>0.75:
            text = "This sequence has a GC content higher than 75%.<br>In the data of Tsai et al Nat Biotech 2015, the two guide sequences with a high GC content had almost as many off-targets as all other sequences combined. We do not recommend using guide sequences with such a high GC content."
            print "<br>"
            htmlWarn(text)
            print ' High GC content<br>'

        if gcContent(guideSeq)<0.25:
            text = "This sequence has a GC content lower than 25%.<br>In the data of Wang/Sabatini/Lander Science 2014, guides with a very low GC content had low cleavage efficiency."
            print "<br>"
            htmlWarn(text)
            print ' Low GC content<br>'

        print "<br>"
        if len(mutEnzymes)!=0:
            print "Restr. Enzymes:"
            print ", ".join([x for x,y,z in mutEnzymes.keys()])
            htmlHelp("Restriction enzymes potentially useful for screening mutations induced by the guide RNA.<br>These enzyme sites overlap the main cleavage site 3bp 5' to the PAM.<br>Digestion of the PCR product with this enzyme will not cut the product if the genome was mutated by Cas9. This is a lot easier than screening with the T7 assay, Surveyor or sequencing.")
        print "</small>"
        print "</td>"

        # off-target score, aka specificity score aka MIT score
        if not cpf1Mode:
            print "<td>"
            if guideScore==None:
                print "No matches"
            else:
                print "%d" % guideScore
            print "</td>"

        # eff scores
        if not cpf1Mode:
            if effScores==None:
                print '<td colspan="%d">Too close to end</td>' % len(scoreNames)
                htmlHelp("The efficiency scores require some flanking sequence<br>This guide does not have enough flanking sequence in your input sequence and could not be extended as it was not found in the genome.<br>")
            else:
                for scoreName in scoreNames:
                    # out-of-frame and prox. gc need special treatment
                    if scoreName in ["oof", "proxGc"]:
                        continue
                    score = effScores.get(scoreName, None)
                    if score!=None:
                        effScoresCount += 1
                    if score==None:
                        print '''<td>--</td>'''
                    elif scoreName=="ssc":
                        # save some space
                        numStr = '%.1f' % (float(score))
                        print '''<td>%s</td>'''  % numStr
                    elif scoreDigits.get(scoreName, 0)==0:
                        print '''<td>%d</td>''' % int(score)
                    else:
                        print '''<td>%0.1f</td>''' % (float(score))
                #print "<!-- %s -->" % seq30Mer

            if showProxGcCol:
                print "<td>"
                # close GC > 4
                finalGc = int(effScores.get("finalGc6", -1))
                if finalGc==1:
                    print "+"
                elif finalGc==0:
                    print "-"
                else:
                    print "--"

                # main motif is "NGG" and last nucleotides are GGNGG
                #if pam=="NGG" and patMatch(guideSeq[-2:], "GG"):
                if int(effScores.get("finalGg", 0))==1:
                    print "<br>"
                    print "<small>-GG</small>"
                print "</td>"

        print "<td>"
        oofScore = str(effScores.get("oof", None))
        if oofScore==None:
            print "--"
        else:
            print """<a href="%s?batchId=%s&pamId=%s&showMh=1" target=_blank class="tooltipster" title="The score indicates how likly out-of-frame deletions are. Click to show the induced deletions based on the micro-homology around the cleavage site.">%s</a>""" % (myName, batchId, pamId, oofScore)
            #print """<br><br><small><a href="%s?batchId=%s&pamId=%s&showMh=1" target=_blank class="tooltipster">Micro-homology</a></small>""" % (myName, batchId, pamId)
        print "</td>"

        # mismatch description
        print "<td>"
        #otCount = sum([int(x) for x in otDesc.split("/")])
        if otData==None:
            # no genome match
            print otDesc
            htmlHelp("Sequence was not found in genome.<br>If you have pasted a cDNA sequence, note that sequences that overlap a splice site cannot be used as guide sequences<br>This warning also appears if you have selected the wrong or no genome.")
        elif subOptMatchCount > MAXOCC:
            print ("Repeat")
            htmlHelp("At <= 4 mismatches, %d hits were found in the genome for this sequence. <br>This guide is a repeated region, it is too unspecific.<br>Usually, CRISPR cannot be used to target repeats." % subOptMatchCount)
        else:
            print otDesc
            print "<br>"

            # mismatch description, last 12 bp
            print '<small style="color:grey">'+last12Desc+"</small><br>"
            otCount = len(otData)
            print "<br><small>%d off-targets</small>" % otCount
        print "</td>"

        # links to offtargets
        print "<td><small>"
        if otData!=None:
            if len(otData)>500 and len(guideData)>1:
                otData, cutoff = findOtCutoff(otData)
                if cutoff!=None:
                    print "Too many off-targets. Showing those with score &gt;%0.1f " % cutoff
                else:
                    print "Too man off-targets. Showing "+str(len(otData))

                htmlHelp("This guide sequence has a high number of off-targets, its use is discouraged.<br>To show all off-targets, paste only the guide sequence into the input sequence box.")

            otLinks = makeOtBrowserLinks(otData, chrom, dbInfo, pamId)

            print "\n".join(otLinks[:3])
            if len(otLinks)>3:
                cssPamId = pamId.replace("-","minus").replace("+","plus") # +/-: not valid in css
                cssPamId = cssPamId+"More"
                print '<div id="%s" class="otMore" style="display:none; width:100%%">' % cssPamId
                print "\n".join(otLinks[3:])

                print '''<a style="float:right;text-decoration:underline" href="%s?batchId=%s&pamId=%s&otPrimers=1" id="%s">''' % (myName, batchId, pamId, cssPamId)
                print 'Off-target primers</a>'

                print '</div>'

                print '''<a id="%sMoreLink" class="otMoreLink" onclick="showAllOts('%s')">''' % (cssPamId, cssPamId)
                print 'show all...</a>'


                print '''<a id="%sLessLink" class="otLessLink" style="display:none" onclick="showLessOts('%s')">''' % (cssPamId, cssPamId)
                print 'show less...</a>'



        print "</small></td>"

        print "</tr>"
        count = count+1

    print "</table>"
    printDownloadTableLinks(batchId)

    printNoEffScoreFoundWarn(effScoresCount)

def linkLocalFiles(listFname):
    """ write a <link> statement for each filename in listFname. Version them via mtime
    (-> browser cache)
    """
    for fname in open(listFname).read().splitlines():
        fname = fname.strip()
        if not isfile(fname):
            fname = join(HTMLDIR, fname)
            if not isfile(fname):
                print "missing: %s<br>" % fname
                continue
        mTime = str(os.path.getmtime(fname)).split(".")[0] # seconds is enough
        if fname.endswith(".css"):
            #url = fname.replace("/var/www/", "http://tefor.net/")
            print "<link rel='stylesheet' media='screen' type='text/css' href='%s%s?%s'/>" % (HTMLPREFIX, fname, mTime)

def printHeader(batchId, title):
    " print the html header "

    print "<html><head>"

    if title==None:
        if batchName!="":
            print """<title>CRISPOR - %s</title>""" % batchName
        else:
            print """<title>CRISPOR</title>"""
    else:
        print """<title>%s</title>""" % title
    print """
<meta name='description' content='Design CRISPR guides with off-target and efficiency predictions, for more than 100 genomes.'/>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />
<meta property='fb:admins' content='692090743' />
<meta name="google-site-verification" content="OV5GRHyp-xVaCc76rbCuFj-CIizy2Es0K3nN9FbIBig" />
<meta property='og:type' content='website' />
<meta property='og:url' content='http://tefor.net/crispor/crispor.py' />
<meta property='og:image' content='http://tefor.net/crispor/image/CRISPOR.png' />

"""

    # load jquery from local copy, not from CDN, for offline use
    print("""<script src='%sjs/jquery.min.js'></script>
<script src='%sjs/jquery-ui.min.js'></script>
""" % (HTMLPREFIX, HTMLPREFIX))

    linkLocalFiles("includes.txt")

    print '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster-shadow.css" />' % HTMLPREFIX

    # the UFD combobox, https://code.google.com/p/ufd/wiki/Usage
    # patched to allow mouse wheel
    # https://code.google.com/p/ufd/issues/detail?id=86&q=mouse%20wheel
    print '<script type="text/javascript" src="%sjs/jquery.ui.ufd.js"></script>' % HTMLPREFIX
    #print '<link rel="stylesheet" type="text/css" href="%sstyle/ufd-base.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css" href="%sstyle/plain.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css"  href="%sstyle/jquery-ui.css" />' % HTMLPREFIX
    print '<script type="text/javascript" src="js/jquery.tooltipster.min.js"></script>'

    # override the main TEFOR css
    print '<style>'
    print 'body { text-align: left; float: left} '
    print '.contentcentral { text-align: left; float: left} '
    print '</style>'

    # activate tooltipster
   #theme: 'tooltipster-shadow',
    print ("""
    <script> 
    $(document).ready(function() { 
        $('.tooltipster').tooltipster({ 
            minWidth: 0,
            contentAsHTML: true,
            maxWidth:400,
            arrow: false,
            interactive: true,
            speed : 0
        }); });
    $(document).ready(function() {
        $('.tooltipsterInteract').tooltipster({
            minWidth: 0,
            contentAsHTML: true,
            maxWidth:400,
            interactive: true,
            onlyOne: true,
            arrow: false,
            speed : 0
        }); });
    </script> """)

    # activate jqueryUI tooltips
    print ("""
    <script>
    $(function () {
       $(".tooltip").tooltip({
       relative : true,
       tooltipClass : "alignStyle",
       content: function () {
       return '<div style="width:300px">'+$(this).prop('title')+"</div>";
       }
      });
    });

    $(function () {
       $(".tooltipAuto").tooltip({
       contentAsHtml : true
      });
    });
    </script>""")


    # style of Jquery UI tooltips, default style is div.ui-tooltip
    print("""<style>
        .alignStyle {
            background-color: #FFFFFF;
            width: 350px;
            max-width: 400px;
            height: 110px;
            position : absolute;
            text-align: left;
            border:1px solid #cccccc;
        }
            </style>""")

    # style from https://css-tricks.com/rotated-table-column-headers/ to rotate table headers
    print("""<style>
       th.rotate {
         /* Something you can count on */
         /* height: 10px; */
         white-space: nowrap;
       }
       
       th.rotate > div {
         float:left;
         white-space: nowrap;
         position: relative;
         border-style: none;
      """)

    # if we're showing all scores, we have very little space, so turn by 
    # 90degrees. otherwise, we can afford 45 degrees, which is easier to read.
    #if cgiParams.get("showAllScores", "0")=="1":
    print("""
    -webkit-transform: rotate(-90);
    -moz-transform: rotate(270deg);
    -ms-transform: rotate(270deg);
    -o-transform: rotate(270deg);
    transform: rotate(270deg);
    width: 25px;""")

    #else:
          #print("""
         #-webkit-transform: rotate(-45deg);
         #-moz-transform: rotate(315deg);
         #-ms-transform: rotate(315deg);
         #-o-transform: rotate(315deg);
         #transform: rotate(315deg);
         #width: 25px;""")

    print("""
       }

       th.rotate > div > span {
         /* border-bottom: 1px solid #ccc; */
         padding: 0px 3px;
         white-space: nowrap;
       }
    </style>""")


    print("</head>")

    print'<body id="wrapper">'
    
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
    for pamId, pamStart, guideStart, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, len(pam), True):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()

def runCmd(cmd, ignoreExitCode=False):
    " run shell command, check ret code, replaces BIN and SCRIPTS special variables "
    cmd = cmd.replace("$BIN", binDir)
    cmd = cmd.replace("$SCRIPT", scriptDir)
    cmd = "set -o pipefail; " + cmd
    debug("Running %s" % cmd)
    ret = subprocess.call(cmd, shell=True, executable="/bin/bash")
    if ret!=0 and not ignoreExitCode:
        if commandLineMode:
            logging.error("Error: could not run command %s." % cmd)
            sys.exit(1)
        else:
            print "Server error: could not run command %s, error %d.<p>" % (cmd, ret)
            print "please send us an email, we will fix this error as quickly as possible. %s " % contactEmail
            sys.exit(0)


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
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    startDict, endSet = findAllPams(seq, pam)
    pamInfo = list(flankSeqIter(seq, startDict, len(pam), False))

    guideIds = []
    guides = []
    longSeqs = []
    for pamId, startPos, guideStart, strand, guideSeq, pamSeq in pamInfo:
        guides.append(guideSeq+pamSeq)
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam))
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

def findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname):
    " align faFname to genome and create matchedBedFname "
    matchesBedFname = batchBase+".matches.bed"
    saFname = batchBase+".sa"
    pamLen = len(pam)
    genomeDir = genomesDir # make var local, see below

    open(matchesBedFname, "w") # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode
    #if len(parseFasta(open(faFname)))==1 and not commandLineMode:
        #global MAXOCC
        #global maxMMs
        #MAXOCC=max(HIGH_MAXOCC, MAXOCC)
        #maxMMs=HIGH_maxMMs

    maxDiff = maxMMs
    queue.startStep(batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff)
    convertMsg = "Converting alignments"
    seqLen = GUIDELEN

    bwaM = MFAC*MAXOCC # -m is queue size in bwa
    cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(saFname)s" % locals()
    runCmd(cmd)

    queue.startStep(batchId, "saiToBed", convertMsg)
    maxOcc = MAXOCC # create local var from global
    # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
    # the sorting should improve the twoBitToFa runtime
    cmd = "$BIN/bwa samse -n %(maxOcc)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | $SCRIPT/samToBed %(pam)s | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomeDir)s/%(genome)s/%(genome)s.sizes stdout >> %(matchesBedFname)s " % locals()
    runCmd(cmd)

    # arguments: guideSeq, mainPat, altPats, altScore, passX1Score
    filtMatchesBedFname = batchBase+".filtMatches.bed"
    queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
    altPats = ",".join(offtargetPams.get(pam, ["na"]))
    bedFnameTmp = bedFname+".tmp"
    altPamMinScore = str(ALTPAMMINSCORE)
    # EXTRACTION OF SEQUENCES + ANNOTATION
    # twoBitToFa was 15x slower than python's twobitreader, after markd's fix it should be OK
    cmd = "$BIN/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -bed=%(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
    #cmd = "$SCRIPT/twoBitToFaPython %(genomeDir)s/%(genome)s/%(genome)s.2bit %(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
    runCmd(cmd)

    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # if we have gene model segments, annotate them, otherwise just use the chrom position
    if isfile(segFname):
        queue.startStep(batchId, "genes", "Annotating matches with genes")
        cmd = "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 2> %(batchBase)s.log > %(bedFnameTmp)s " % locals()
        runCmd(cmd)
    else:
        queue.startStep(batchId, "chromPos", "Annotating matches with chromosome position")
        annotateBedWithPos(filtMatchesBedFname, bedFnameTmp)

    # make sure the final bed file is never in a half-written state, 
    # as it is our signal that the job is complete
    shutil.move(bedFnameTmp, bedFname)
    queue.startStep(batchId, "done", "Job completed")
    return bedFname

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

def findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pamPat, bedFname):
    " align guides with pam in faFname to genome and write off-targets to bedFname "
    tmpDir = batchBase+".bowtie.tmp"
    os.mkdir(tmpDir)

    # make sure this directory gets removed, no matter what
    global tmpDirsDelExit
    tmpDirsDelExit.append(tmpDir)
    if not DEBUG:
        atexit.register(delTmpDirs)

    # write out the sequences for bowtie
    queue.startStep(batchId, "seqPrep", "preparing sequences")
    bwFaFname = abspath(join(tmpDir, "bowtieIn.fa"))
    guideSeqs, qSeqs, allPamSeqs = writeBowtieSequences(faFname, bwFaFname, pamPat)

    genomePath =  abspath(join(genomesDir, genome, genome))
    oldCwd = os.getcwd()

    # run bowtie
    queue.startStep(batchId, "bowtie", "aligning with bowtie")
    os.chdir(tmpDir) # bowtie writes to hardcoded output filenames with --refout
    # -v 3 = up to three mismatches
    # -y   = try hard
    # -t   = print time it took
    # -k   = output up to X alignments
    # -m   = do not output any hit if a read has more than X hits
    # --max = write all reads that exceed -m to this file
    # --refout = output in bowtie format, not SAM
    # --maxbts=2000 maximum number of backtracks
    # -p 4 = use four threads
    # --mm = use mmap
    maxOcc = MAXOCC # meaning in BWA: includes any PAM, in bowtie we have the PAM in the input sequence
    cmd = "$BIN/bowtie %(genomePath)s -f %(bwFaFname)s  -v 3 -y -t -k %(maxOcc)d -m %(maxOcc)d dummy --max tooManyHits.txt --mm --refout --maxbts=2000 -p 4" % locals()
    runCmd(cmd)
    os.chdir(oldCwd)

    queue.startStep(batchId, "parse", "parsing alignments")
    pamLen = len(pamPat)
    hits = parseRefout(tmpDir, guideSeqs, pamLen)

    queue.startStep(batchId, "scoreOts", "scoring off-targets")
    # make the list of alternative PAM sequences
    altPats = offtargetPams.get(pamPat, [])
    altPamSeqs = []
    for altPat in altPats:
        altPamSeqs.extend(expandIupac(altPat))

    # iterate over bowtie hits and write to a BED file with scores
    # if the hit looks OK (right PAM + score is high enough)
    tempBedPath = join(tmpDir, "bowtieHits.bed")
    tempFh = open(tempBedPath, "w")

    offTargets = {}
    for guideIdWithMod, chrom, start, end, strand, _, tSeq in hits:
        guideId = guideIdWithMod.split(".")[0]
        guideSeq = guideSeqs[guideId]
        genomePamSeq = tSeq[-pamLen:]
        logging.debug( "PAM seq: %s of %s" % (genomePamSeq, tSeq))
        if genomePamSeq in altPamSeqs:
            minScore = ALTPAMMINSCORE
        elif genomePamSeq in allPamSeqs:
            minScore = MINSCORE
        else:
            logging.debug("Skipping off-target for %s: %s:%d-%d" % (guideId, chrom, start, end))
            continue

        logging.debug("off-target minScore = %f" % minScore )

        # check if this match passes the off-target score limit
        if pamPat=="TTTN":
            otScore = 0.0
        else:
            tSeqNoPam = tSeq[:-pamLen]

            otScore = calcHitScore(guideSeq, tSeqNoPam)

            if otScore < minScore:
                logging.debug("off-target not accepted")
                continue

        editDist = getEditDist(guideSeq, tSeqNoPam)
        guideHitCount = 0
        guideId = guideId.split(".")[0] # full guide ID looks like s33+.0:A>T
        name = guideId+"|"+strand+"|"+str(editDist)+"|"+tSeq+"|"+str(guideHitCount)+"|"+str(otScore)
        row = [chrom, str(start), str(end), name]
        # this way of collecting the features will remove the duplicates
        otKey = (chrom, start, end, strand, guideId)
        logging.debug("off-target key is %s" % str(otKey))
        offTargets[ otKey ] = row

    for rowKey, row in offTargets.iteritems():
        tempFh.write("\t".join(row))
        tempFh.write("\n")

    tempFh.flush()

    # create a tempfile which is moved over upon success
    # makes sure we do not leave behind a half-written file if 
    # we crash later
    tmpFd, tmpAnnotOffsPath = tempfile.mkstemp(dir=tmpDir, prefix="annotOfftargets")
    tmpFh = open(tmpAnnotOffsPath, "w")

    # get name of file with genome locus names
    genomeDir = genomesDir # make local var
    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # annotate with genome locus names
    cmd = "$BIN/overlapSelect %(segFname)s %(tempBedPath)s stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(tmpAnnotOffsPath)s" % locals()
    runCmd(cmd)

    shutil.move(tmpAnnotOffsPath, bedFname)
    queue.startStep(batchId, "done", "Job completed")

    if DEBUG:
        logging.info("debug mode: Not deleting %s" % tmpDir)
    else:
        shutil.rmtree(tmpDir)

def processSubmission(faFname, genome, pam, bedFname, batchBase, batchId, queue):
    """ search fasta file against genome, filter for pam matches and write to bedFName 
    optionally write status updates to work queue.
    """
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
        findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pam, bedFname)
    else:
        findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname)

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
    global allGenomes
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

def printOrgDropDown(lastorg, genomes):
    " print the organism drop down box. "
    print '<select id="genomeDropDown" class style="max-width:400px" name="org" tabindex="2">'
    print '<option '
    if lastorg == "noGenome":
        print 'selected '
    print 'value="noGenome">-- No Genome: no specificity, only cleavage efficiency scores (max. len 25kbp)</option>'

    for db, desc in genomes:
        print '<option '
        if db == lastorg :
            print 'selected '
        print 'value="%s">%s</option>' % (db, desc)

    print "</select>"
    #print ('''
      #<script type="text/javascript">
      #$("#genomeDropDown").ufd({maxWidth:350, listWidthFixed:false});
      #</script>''')
    print ('''<br>''')

def printPamDropDown(lastpam):
    
    print '<select style="float:left" name="pam" tabindex="3">'
    for key,value in pamDesc:
        print '<option '
        if key == lastpam :
            print 'selected '
        print 'value="%s">%s</option>' % (key, value)
    print "</select>"

def printForm(params):
    " print html input form "
    scriptName = basename(__file__)

    genomes = readGenomes()

    haveHuman = False
    for g in genomes:
        if g[0]=="hg19":
             haveHuman = True

    # The returned cookie is available in the os.environ dictionary
    cookies=Cookie.SimpleCookie(os.environ.get('HTTP_COOKIE'))
    if "lastorg" in cookies and "lastseq" in cookies and "lastpam" in cookies:
       lastorg   = cookies['lastorg'].value
       lastseq   = cookies['lastseq'].value
       lastpam   = cookies['lastpam'].value
    else:
       if not haveHuman:
           global DEFAULTSEQ
           global DEFAULTORG
           DEFAULTSEQ = ALTSEQ
           DEFAULTORG = ALTORG
       lastorg = DEFAULTORG
       lastseq = DEFAULTSEQ
       lastpam = DEFAULTPAM

    print """
<form id="main-form" method="post" action="%s">

<!-- <strong>Web server maintenance at TEFOR.NET: Jan 13 - Jan 15 2016<br>
Site temporarily moved to UCSC. The performance here is somewhat slower.<br>
Site should be back online at the original URL during Jan 16 2016<p></strong> -->

 <div style="text-align:left; margin-left: 50px">
 CRISPOR is a program that helps design, evaluate and clone guide sequences for the CRISPR/Cas9 system.

<span class="introtext">
    <div onclick="$('.about-us').toggle('fast');" class="title" style="cursor:pointer;display:inline;font-size:large;font-style:normal">
        <img src="%simage/info-small.png" style="vertical-align:text-top;">
    </div>
    <div class="about-us"><br>
    CRISPOR v3.0 uses the BWA algorithm to identify guide RNA sequences for CRISPR mediated genome editing.<br>
    It searches for off-target sites (with and without mismatches), shows them in a table and annotates them with flanking genes.<br>
    For more information on principles of CRISPR-mediated genome editing, check the <a href="https://www.addgene.org/CRISPR/guide/">Addgene CRISPR guide</a>.</div>
</span>

<br><i>Wondering which score is best for you? Have a look at the <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2">CRISPOR paper</a><small> (and email us any questions)</small></i>

 </div>

<div class="windowstep subpanel" style="width:50%%;">
    <div class="substep">
        <div class="title">
            Step 1
        </div>
            
        Sequence name (optional): <input type="text" name="name" size="20"><br>

        Enter a single genomic sequence, &lt; %d bp, typically an exon
        <img src="%simage/info-small.png" title="CRISPOR conserves the lowercase and uppercase format of your sequence, allowing to highlight sequence features of interest such as ATG or STOP codons.<br>Avoid using cDNA sequences as input, CRISPR guides that straddle splice sites are unlikely to work.<br>You can paste a single 20bp guide and even multiple guides, separated by N characters." class="tooltipster">
    <br>
    <small><a href="javascript:clearInput()">Clear Box</a> - </small>
    <small><a href="javascript:resetToExample()">Reset to default</a></small>
    </div>

    <textarea tabindex="1" style="width:100%%" name="seq" rows="12"
              placeholder="Paste here the genomic - not cDNA - sequence of the exon you want to target. Maximum size %d bp.">%s</textarea>
      <small>Text case is preserved, e.g. you can mark ATGs with lowercase.<br>Instead of a sequence, you can paste a chromosome range, e.g. chr1:11908-12378</small>
</div>
<div class="windowstep subpanel" style="width:40%%">
    <div class="substep" style="margin-bottom: 1px">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep2').toggle('fast')">
            Step 2
        </div>
        Select a genome
    </div>
    """% (scriptName, HTMLPREFIX, MAXSEQLEN, HTMLPREFIX, MAXSEQLEN, lastseq)

    printOrgDropDown(lastorg, genomes)
    print """
    <div id="trackHubNote" style="margin-bottom:5px">
    <small>Note: we have pre-calculated all exonic guides for this species. Have a look at the <a id='hgTracksLink' target=_blank href="">UCSC Genome Browser</a>.</small>
    </div>
    """
    print '<small style="float:left">Missing a genome? Send us <a href="mailto:%s">email</a></small>' % (contactEmail)
    print """
    </div>
    <div class="windowstep subpanel" style="width:40%%; height:158px">
    <div class="substep">
    <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
        Step 3
        <img src="%simage/info-small.png" title="The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href='http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html' target='_blank'>Kleinstiver et al</a>, Nature 2015" class="tooltipsterInteract">
        </div>
        Select a Protospacer Adjacent Motif (PAM)
    </div>
    """ % HTMLPREFIX

    printPamDropDown(lastpam)
    print """
    <div style="width:40%; margin-top: 50px; margin-left:50px; text-align:center; display:block">
    <input type="submit" name="submit" value="SUBMIT" tabindex="4"/>
    </div>
    </div>
    """
    print """
<script>
/* set the dropbox to hg19 and paste the example sequence into the input box. */
function resetToExample() {
    $("textarea[name='seq']").val("%s");
    $("#genomeDropDown").val("%s");
    $("select[name='pam']").val("NGG");
    }

/* clear the sequence input box */
function clearInput() {
    $("textarea[name='seq']").val("");
    }

</script>
<script>
    /* hide the track hub note if genome is not hg19 */
    ucscTrackDbs=['hg19', 'hg38', 'rn5', 'mm10', 'mm9', 'ci2', 'danRer7', 'sacCer3', 'dm6'];
    function showHideHubNote() {
        var valSel = $("#genomeDropDown").val();
        if (jQuery.inArray(valSel, ucscTrackDbs)!=-1)
            {
            $("#trackHubNote").css('visibility', 'visible');
            $("#hgTracksLink").attr("href", "http://genome-test.soe.ucsc.edu/cgi-bin/hgTracks?db="+valSel+"&crispr=show");
            }
        else
            $("#trackHubNote").css('visibility', 'hidden');
    }
    $("#genomeDropDown").on('change', showHideHubNote);
    showHideHubNote();
</script>

</form>
    """ % (DEFAULTSEQ, DEFAULTORG)

def readBatchParams(batchId):
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

def newBatch(batchName, seq, org, pam, skipAlign=False):
    """ obtain a batch ID and write seq/org/pam to their files.
    Return batchId, position string and a 100bp-extended seq, if possible.
    """
    batchId = makeTempBase(seq, org, pam, batchName)
    if skipAlign:
        chrom, start, end, strand = None, None, None, None
    else:
        chrom, start, end, strand = findBestMatch(org, seq)
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
        extSeq = extendAndGetSeq(org, chrom, start, end, strand)
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

def printQueryNotFoundNote(dbInfo):
    print "<div class='title'>Query sequence, not present in the genome of %s</div>" % dbInfo.scientificName
    print "<div class='substep'>"
    print "<em><strong>Note:</strong>The query sequence was not found in the selected genome."
    print "This can be a valid query, e.g. a GFP sequence.<br>"
    print "If not, you might want to check if you selected the right genome for your query sequence.<br>"
    print "When reading the list of guide sequences and off-targets below, bear in mind that the software cannot distinguish off-targets from on-targets now, so some 0-mismatch targets are expected. In this case, the scores of guide sequences are too low.<br>"
    print "Because there is no flanking sequence available, the guides in your sequence that are within 50bp of the ends will have no efficiency scores. The efficiency scores will instead be shown as '--'. Include more flanking sequence > 50bp to obtain the scores.<p>"
    print "</em></div>"

def getOfftargets(seq, org, pam, batchId, startDict, queue):
    """ write guides to fasta and run bwa or use cached results.
    Return name of the BED file with the matches.
    Write progress status updates to queue object.
    """
    batchBase = join(batchDir, batchId)
    otBedFname = batchBase+".bed"
    flagFile = batchBase+".running"

    if isfile(flagFile):
       errAbort("This sequence is still being processed. Please wait for ~20 seconds "
           "and try again, e.g. by reloading this page. If you see this message for "
           "more than 2-3 minutes, please send an email %s.net. Thanks!" % contactEmail)

    if not isfile(otBedFname) or commandLineMode:
        # write potential PAM sites to file 
        faFname = batchBase+".fa"
        writePamFlank(seq, startDict, pam, faFname)
        if commandLineMode:
            processSubmission(faFname, org, pam, otBedFname, batchBase, batchId, queue)
        else:
            # umask is not respected by sqlite, bug http://www.mail-archive.com/sqlite-users@sqlite.org/msg59080.html
            q = JobQueue()
            ip = os.environ["REMOTE_ADDR"]
            wasOk = q.addJob("search", batchId, "ip=%s,org=%s,pam=%s" % (ip, org, pam))
            if not wasOk:
                #print "CRISPOR job is running..." % batchId
                pass
            q.close()
            return None

    return otBedFname

def startAjaxWait(batchId):
    """ print the ajax script to stdout """
    scriptName = basename(__file__)
    #$.ajax({ url: "server", success: function(data){
    #//$.getJSON( "%(scriptName)s?batchId=%(batchId)s&ajaxStatus=1", gotStatus);
    print """
    Status: <div id="statusEl">Waiting</div>
    <script>
    function gotStatus( data ) {
      status = data["status"];
      if (status==null || status=="null")
          { window.location.href=window.location.href; }
      else
          { $("#statusEl").text(status); }
      };
     
    (function poll() {
       setTimeout(function(){
            $.getJSON("%(scriptName)s?batchId=%(batchId)s&ajaxStatus=1", gotStatus);
            poll();
          } , 10000)})();

    </script>
    """ % locals()

def showPamWarning(pam):
    if pam=="TTTN":
        print '<div style="text-align:left">'
        print "<strong>Note:</strong> You are using the Cpf1 enzyme."
        print "Note that currently there are no on- or off-target scores for Cpf1 so no scores are shown below and the guides are not sorted."
        print '</div>'
    elif pam!="NGG":
        print '<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">'
        print "<strong>Note:</strong> Your query involves a Cas9 that is not from S. Pyogenes. "
        print "Please bear in mind that specificity and efficiency scores were designed using data with S. Pyogenes Cas9 and might not be applicable to this particular Cas9.<br>"
        print '</div>'

def showNoGenomeWarning(dbInfo):
    if dbInfo==None:
        print('<div style="text-align:left;"><strong>Note:</strong> As there is no genome that can be used to get flanking sequence for your sequence, efficiency scores 50bp from the start or the end of your sequence cannot be calculated and are shown as "--". If needed, extend the input sequence and retry.</div>')

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

def printStatus(batchId):
    " print status, not using any Ajax "
    q = JobQueue()
    status = q.getStatus(batchId)
    q.close()

    errorState = False

    if "Traceback" in status:
        print "<!--"
        print status
        print "-->"
        status = "An error occured during processing.<br> Please send an email to %s and tell us that the failing batchId was %s.<br>We can usually fix this quickly. Thanks!" % (contactEmail, batchId)
        errorState = True
    else:
        print('<meta http-equiv="refresh" content="10" >')
        print("CRISPOR job has been submitted.<p>")

    if status==None:
        status = "Batch completed. Refresh page to show results."

    print("Job Status: <tt>%s</tt><p>" % status)

    if not errorState:
        print("<small>This page will refresh every 10 seconds</small>")

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
            #afList = infoDict["AF"].split(",")
            #altAllList = altAll.split(",")
            #newAltAllList = []
            #newAfList = []
            #for af, altAll in zip(afList, altAllList):
                #afNum = float(af)
                #if afNum < minFreq:
                    #continue
                #newAltAllList.append(altAll)
                #newAfList.append(af)
            ##if len(newAltAllList)==0:
                #continue
            #altAll = ",".join(newAltAllList)
            #infoDict["AF"] = ",".join(newAfList)
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
        #infoDict["dbg"] = "%s:rel=%d:abs=%s:%s->%s" % (chrom, relPos, varPos, refAll, altAll)
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
        print ("<script>")
        print ('''history.replaceState('crispor.py', document.title, '?batchId=%s');''' % (batchId))
        print ("</script>")

    setupPamInfo(pam)

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

    if otBedFname is None:
        # this can happen only in CGI mode. Job has been added to the queue or is not done yet. 
        #startAjaxWait(batchId)
        printStatus(batchId)
        return

    # be more sensitive if only a single guide seq is run
    #if len(startDict)==1:
        #global MAXOCC
        #MAXOCC=max(HIGH_MAXOCC, MAXOCC)

    if dbInfo==None:
        print "<div class='title'>No Genome selected, specificity scoring is deactivated</div>"
        print('<div style="text-align:left;"><strong>Note:</strong> There are no predicted off-targets below and all specificity scores are shown in red as their score is 0. <br></div>')

    elif position=='?':
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
        print makeBrowserLink(dbInfo, genomePosStr, oneBasedPosition, mouseOver, ["tooltipster"])+"</span>, "
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

    varHtmls = varDictToHtml(varDict, seq)
    showSeqAndPams(seq, startDict, pam, guideScores, varHtmls, varDbs, varDb, minFreq)

    showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, chrom, varHtmls)

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

def findBestMatch(genome, seq):
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
    runCmd(cmd)

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
    pamInfo = list(flankSeqIter(inSeq, startDict, len(pam), False))
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

    global batchName
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
    
def parseArgs():
    " parse command line options into args and options "
    parser = optparse.OptionParser("""usage: %prog [options] org fastaInFile guideOutFile 

Command line interface for the Crispor tool.

    org          = genome identifier, like hg19 or ensHumSap
    fastaInFile  = Fasta file
    guideOutFile = tab-sep file, one row per guide

    If many guides have to scored in batch: Add GGG to them to make them valid
    guides, separate these sequences by N characters and supply as a single
    fasta sequence, a few dozen to ~100 per file.
    """)

    parser.add_option("-d", "--debug", dest="debug", \
        action="store_true", help="show debug messages, do not delete temp directory")
    parser.add_option("-t", "--test", dest="test", \
        action="store_true", help="run internal tests")
    pamNames = (",".join([x for x,y in pamDesc]))
    parser.add_option("-p", "--pam", dest="pam", \
        action="store", help="PAM-motif to use, default %default. TTTN triggers special Cpf1 behavior: no scores anymore + the PAM is assumed to be 5' of the guide. Common PAMs are: " + pamNames, default="NGG")
    parser.add_option("-o", "--offtargets", dest="offtargetFname", \
        action="store", help="write offtarget info to this filename")
    parser.add_option("-m", "--maxOcc", dest="maxOcc", \
        action="store", type="int", help="MAXOCC parameter, guides with more matches are excluded")
    parser.add_option("", "--mm", dest="mismatches", \
        action="store", type="int", help="maximum number of mismatches, default %default", default=4)
    parser.add_option("", "--bowtie", dest="bowtie", \
        action="store_true", help="new: use bowtie as the aligner. Do not use. Bowtie misses many off-targets.")
    parser.add_option("", "--skipAlign", dest="skipAlign", \
        action="store_true", help="do not align the input sequence. The on-target will be a random match with 0 mismatches.")
    parser.add_option("", "--noEffScores", dest="noEffScores", \
        action="store_true", help="do not calculate the efficiency scores")
    parser.add_option("", "--minAltPamScore", dest="minAltPamScore", \
        action="store", type="float", help="minimum MIT off-target score for alternative PAMs, default %default", \
        default=ALTPAMMINSCORE)
    parser.add_option("", "--worker", dest="worker", \
        action="store_true", help="Run as worker process: watches job queue and runs jobs")
    #parser.add_option("", "--user", dest="user", \
        #action="store", help="for the --worker option: switch to this user at program start")
    parser.add_option("", "--clear", dest="clear", \
        action="store_true", help="clear the worker job table and exit")
    parser.add_option("-g", "--genomeDir", dest="genomeDir", \
        action="store", help="directory with genomes, default %default", default=genomesDir)

    #parser.add_option("-f", "--file", dest="file", action="store", help="run on file") 
    (options, args) = parser.parse_args()

    if len(args)==0 and not options.test and not options.worker and not options.clear:
        parser.print_help()
        sys.exit(0)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options

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
    global tmpDirsDelExit
    logging.debug("Removing tmpDirs: %s" % ",".join(tmpDirsDelExit))
    for tmpDir in tmpDirsDelExit:
        if isdir(tmpDir):
            shutil.rmtree(tmpDir)
    tmpDirsDelExit = []



def handleOptions(options):
    " set glpbal vars based on options "
    if options.test:
        runTests()
        import doctest
        doctest.testmod()
        sys.exit(0)

    if options.debug:
        global DEBUG
        DEBUG = True

    if options.noEffScores:
        global doEffScoring
        doEffScoring = False

    # handle the alignment/filtering options
    if options.maxOcc != None:
        global MAXOCC
        MAXOCC = options.maxOcc
        #HIGH_MAXOCC = options.maxOcc

    if options.minAltPamScore!=None:
        global ALTPAMMINSCORE
        ALTPAMMINSCORE = options.minAltPamScore

    if options.mismatches:
        global maxMMs
        maxMMs = options.mismatches

    if options.bowtie:
        global useBowtie
        useBowtie = True

    if options.pam:
        setupPamInfo(options.pam)

def mainCommandLine():
    " main entry if called from command line "
    global commandLineMode
    commandLineMode = True

    args, options = parseArgs()

    handleOptions(options)
    org, inSeqFname, outGuideFname = args

    skipAlign = False
    if options.skipAlign:
        skipAlign = True

    # different genomes directory?
    if options.genomeDir != None:
        global genomesDir
        genomesDir = options.genomeDir

    # get sequence
    seqs = parseFasta(open(inSeqFname))

    # make a directory for the temp files
    # and put it into a global variable, so all functions will use it
    global batchDir
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
        batchId, position, extSeq = newBatch(seqId, seq, org, pam, skipAlign)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        if position=="?":
            logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pam)
        otBedFname = getOfftargets(seq, org, pam, batchId, startDict, ConsQueue())
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

def sendStatus(batchId):
    " send batch status as json "
    q = JobQueue()
    status = q.getStatus(batchId)
    q.close()

    if status==None:
        d = {"status":status}
    elif "Traceback" in status:
        d = {"status" : "An error occured. Please send an email to %s and tell us that the failing batchId was %s. We can usually fix this quickly. Thanks!" % (contactEmail, batchId)}
    else:
        d = {"status":status}
    print json.dumps(d)

def cleanJobs():
    """ look for flag file cleanJobs in current dir. If present, remove jobs.db. 
    this is the only way to remove the file, as the jobs.db file is owned by apache 
    """
    if isfile("cleanJobs"):
        os.remove(JOBQUEUEDB)
        os.remove("cleanJobs")



def main():
    mainCommandLine()
main()
