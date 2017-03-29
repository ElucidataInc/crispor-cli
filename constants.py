import os,platform

from os.path import join, isfile, basename, dirname, isdir, abspath
# write debug output to stdout
DEBUG = False
#DEBUG = True
# use bowtie for off-target search?
useBowtie = False
# calculate the efficienc scores?
doEffScoring = True
# BWA: allow up to X mismatches
maxMMs=4
# maximum number of occurences in the genome to get flagged as repeats. 
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000
# Highly-sensitive mode (not for CLI mode): 
# MAXOCC is increased in processSubmission() and in the html UI if only one
# guide seq is run
# Also, the number of allowed mismatches is increased to 5 instead of 4
#HIGH_MAXOCC=600000
#HIGH_maxMMs=5
ALTPAMMINSCORE = 1.0

# version of crispor
versionStr = "4.0"

# contact email
contactEmail='crispor@tefor.net'



# system-wide temporary directory
TEMPDIR = os.environ.get("TMPDIR", "/tmp")
#TEMPDIR = "/dev/shm/"



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



DEFAULTORG = 'hg19'
DEFAULTSEQ = 'cttcctttgtccccaatctgggcgcgcgccggcgccccctggcggcctaaggactcggcgcgccggaagtggccagggcgggggcgacctcggctcacagcgcgcccggctattctcgcagctcaccatgGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC'

# used if hg19 is not available
ALTORG = 'sacCer3'
ALTSEQ = 'ATTCTACTTTTCAACAATAATACATAAACatattggcttgtggtagCAACACTATCATGGTATCACTAACGTAAAAGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTTCAGAATATTTCGTACTTACACAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCT'



DEFAULTPAM = 'NGG'


# maximum size of an input sequence 
MAXSEQLEN = 1000
# maximum input size when specifying "no genome"
MAXSEQLEN_NOGENOME = 25000
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
ALTPAMMINSCORE = 1.0
# minimum off-target score of standard off-targets (those that end with NGG)
# This should probably be based on the CFD score these days
# But for now, I'll let the user do the filtering
MINSCORE = 0.0



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
