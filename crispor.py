#!/usr/bin/env python2.7
# the tefor crispr tool
# can be run as a CGI or from the command line

# OOF scores are WRONG for Cpf1! -> where is the cut site?

# python std library
import subprocess, tempfile, logging, atexit, glob, shutil
import Cookie, time, sys, cgi, re, random, platform, os
import hashlib, base64, string, operator, urllib, sqlite3, time
import traceback,  pwd, pickle
import batch_functions
import cleanup_functions
import common_functions
import fasta_functions
import get_offtargets
import parser_functions
import sequence_helper_functions
import scoring_functions
import finishing_functions

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



def debug(msg):
    logging.debug(msg)
    if DEBUG:
        print msg



class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))



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
    if options.bulge_size !=None:
        indel_len = options.bulge_size
        
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
        tmpDirsDelExit=[]
        delBatchDir = cleanup_functions.delBatchDir
        atexit.register(delBatchDir,batchDir,tmpDirsDelExit)

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
        process_parameters = [doEffScoring,useBowtie ,cpf1Mode,batchDir,DEBUG,
                                MAXOCC,ALTPAMMINSCORE,maxMMs,GUIDELEN,addGenePlasmids]
        otBedFname = get_offtargets.getOfftargets(seq, org, pam, batchId, batchDir,startDict, ConsQueue(),indel_len,process_parameters)
        otMatches = get_offtargets.parseOfftargets(otBedFname)
        if options.noEffScores or cpf1Mode:
            effScores = {}
        else:
            effScores = scoring_functions.readEffScores(batchDir,batchId)

        guideData, guideScores, hasNotFound = \
            finishing_functions.mergeGuideInfo(seq, startDict, pam, otMatches, position, effScores,argument_list)

        for row in finishing_functions.iterGuideRows(guideData,addHeaders=True):
            print("row....",row,'\n')
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if options.offtargetFname:
            for row in finishing_functions.iterOfftargetRows(guideData):
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
