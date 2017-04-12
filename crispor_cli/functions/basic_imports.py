import sys

from os.path import join, isfile, dirname




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



# a file crispor.conf in the directory of the script allows to override any global variable
myDir = dirname(__file__)
confPath =join(myDir, "crispor.conf")
if isfile(confPath):
    exec(open(confPath))

cgiParams = None

# ====== END GLOBALS ============