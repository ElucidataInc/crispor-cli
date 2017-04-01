import optparse,logging,sys
import constants as cnst
import test_functions
import common_functions

from os.path import join

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

# minimum off-target score for alternative PAM off-targets
# There is not a lot of data to support this cutoff, but it seems
# reasonable to have at least some cutoff, as otherwise we would show
# NAG and NGA like NGG and the data shows clearly that the alternative
# PAMs are not recognized as well as the main NGG PAM.
# so for now, I just filter out very degenerative ones. the best solution
# would be to have a special penalty on the CFD score, but CFS does not 
# support non-NGG PAMs (is this actually true?)
ALTPAMMINSCORE = 1.0

baseDir = cnst.baseDir
# directory for genomes
genomesDir = join(baseDir, "genomes")

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
    parser.add_option("", "--skipAlign", dest="skipAlign", \
        action="store_true", help="do not align the input sequence. The on-target will be a random match with 0 mismatches.")
    parser.add_option("", "--noEffScores", dest="noEffScores", \
        action="store_true", help="do not calculate the efficiency scores")
    parser.add_option("", "--minAltPamScore", dest="minAltPamScore", \
        action="store", type="float", help="minimum MIT off-target score for alternative PAMs, default %default", \
        default=ALTPAMMINSCORE)
    parser.add_option("-g", "--genomeDir", dest="genomeDir", \
        action="store", help="directory with genomes, default %default", default=genomesDir)

    #parser.add_option("-f", "--file", dest="file", action="store", help="run on file") 
    (options, args) = parser.parse_args()

    if len(args)==0 :
        parser.print_help()
        sys.exit(0)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options

def handleOptions(options):
    " set glpbal vars based on options "
    from constants import DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,useBowtie
    if options.test:
        test_functions.runTests()
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

    # if options.bowtie:
    #     useBowtie = True

    if options.pam:
        GUIDELEN,cpf1Mode,addGenePlasmids = common_functions.setupPamInfo(options.pam)
    argument_list = [DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,
                            useBowtie,GUIDELEN,cpf1Mode,addGenePlasmids]
    
    return argument_list

