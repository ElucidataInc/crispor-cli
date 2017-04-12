import logging,shutil,os,glob
from os.path import isdir,join

def delBatchDir(batchDir,tmpDirsDelExit):
    " called at program exit, for command line mode "
    delTmpDirs(tmpDirsDelExit) # first remove any subdirs
    if not isdir(batchDir):
        return
    logging.debug("Deleting dir %s" % batchDir)
    fnames = glob.glob(join(batchDir, "*"))
    if len(fnames)>50:
        raise Exception("cowardly refusing to remove many temp files")
    for fname in fnames:
        os.remove(fname)
    os.removedirs(batchDir)

def delTmpDirs(tmpDirsDelExit):
    " signal handler at program exit, to remove registered tmp dirs "
    logging.debug("Removing tmpDirs: %s" % ",".join(tmpDirsDelExit))
    for tmpDir in tmpDirsDelExit:
        if isdir(tmpDir):
            shutil.rmtree(tmpDir)
    tmpDirsDelExit = []
