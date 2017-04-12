#!/usr/bin/env python2.7
# the tefor crispr tool
# can be run as a CGI or from the command line

# OOF scores are WRONG for Cpf1! -> where is the cut site?

# python std library
import subprocess, tempfile, logging, atexit, glob, shutil
import Cookie, time, sys, cgi, re, random, platform, os
import hashlib, base64, string, operator, urllib, sqlite3, time
import traceback,  pwd, pickle
from functions import *
import functions.constants as cnst

# from functions.bio_functions import *
from functions.fasta_functions import *
from functions.common_functions import *
from functions.basic_imports import *
from functions.constants import *
from collections import defaultdict, namedtuple
from datetime import datetime
from itertools import product
from os.path import join, isfile,  dirname, isdir, abspath

concatGuideAndPam = sequence_helper_functions.concatGuideAndPam

# ====== FUNCTIONS =====
contentLineDone = False



def debug(msg):
    logging.debug(msg)
    if DEBUG:
        print msg



class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))



def mainCommandLine(args,options):
    " main entry if called from command line "
    # options = parser_functions.parseArgs_package(options)
    opt_keys = options.keys()
    if 'pam' not in opt_keys:
        options['pam']='NGG'
    argument_list = parser_functions.handleOptions(options)
    [DEBUG,doEffScoring,MAXOCC,ALTPAMMINSCORE,maxMMs,
    useBowtie,GUIDELEN,cpf1Mode,addGenePlasmids] = argument_list
    org, inSeqFname, outGuideFname = args
    
    skipAlign = False
    baseDir = cnst.baseDir
    # directory for genomes
    genomesDir = join(baseDir, "genomes")
    indel_len=0

    if 'skipAlign' in opt_keys:
        skipAlign = True
    
    # different genomes directory?
    if 'genomesDir' in opt_keys:
        if options['genomeDir'] != None:
            genomesDir = options['genomeDir']
    if 'bulge_size' in opt_keys:
        if options['bulge_size'] !=None:
            indel_len = options['bulge_size']
        
    # get sequence
    seqs = batch_functions.parseFasta(open(inSeqFname))

    # make a directory for the temp files
    # and put it into a global variable, so all functions will use it
    # global batchDir
    batchDir = tempfile.mkdtemp(dir=TEMPDIR, prefix="crispor")
    logging.debug("Created directory %s" % batchDir)
    if 'debug' in opt_keys:
        logging.info("debug-mode, temporary directory %s will not be deleted" % batchDir)
    else:
        tmpDirsDelExit=[]
        delBatchDir = cleanup_functions.delBatchDir
        atexit.register(delBatchDir,batchDir,tmpDirsDelExit)

    # prepare output files
    guideFh = open(join(batchDir, "guideInfo.tab"), "w")
    guideFh.write("\t".join(guideHeaders)+"\n")
    if 'offtargetFname' in opt_keys:
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
        pam = options['pam']
        batchId, position, extSeq = batch_functions.newBatch(seqId, seq, org, pam,genomesDir,batchDir,skipAlign)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        if position=="?":
            logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pam)
        process_parameters = [doEffScoring,useBowtie ,cpf1Mode,batchDir,DEBUG,
                                MAXOCC,ALTPAMMINSCORE,maxMMs,GUIDELEN,addGenePlasmids]
        otBedFname = get_offtargets.getOfftargets(seq, org, pam, batchId, batchDir,startDict, ConsQueue(),indel_len,process_parameters)
        otMatches = get_offtargets.parseOfftargets(otBedFname)
        if 'noEffScores' in opt_keys or cpf1Mode:
            effScores = {}
        else:
            effScores = scoring_functions.readEffScores(batchDir,batchId)

        guideData, guideScores, hasNotFound = \
       finishing_functions.mergeGuideInfo(seq, startDict, pam, otMatches, position, effScores,argument_list)

        for row in finishing_functions.iterGuideRows(guideData,addHeaders=True):
            print("row....",row,'\n')
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if 'offtargetFname' in opt_keys:
            for row in finishing_functions.iterOfftargetRows(guideData):
                offtargetFh.write("\t".join(row))
                offtargetFh.write("\n")

    guideFh.close()
    shutil.move(guideFh.name, outGuideFname)
    logging.info("guide info written to %s" % outGuideFname)

    if 'offtargetFname' in opt_keys:
        offtargetFh.close()
        shutil.move(offtargetFh.name, options['offtargetFname'])
        logging.info("off-target info written to %s" % options['offtargetFname'])

    if not 'debug'  in opt_keys:
       shutil.rmtree(batchDir)


def main(args,options):
    # options = parser_functions.parseArgs()
    mainCommandLine(args,options)
curdir = os.getcwd()+'/crispor_cli'
args=['sacCer3',curdir+'/input/guide_yeast.fasta',curdir+'/output/yo_guide.tsv']
options = {'offtargetFname':curdir+'/output/yo_off.tsv'}
# main(args,options)