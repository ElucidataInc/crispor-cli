# CRISPOR - a CRISPR/Cas9 assistant 

CRISPOR predicts off-targets in the genome, ranks guides, highlights
problematic guides, designs primers and helps with cloning.  Try it on
http://crispr.org

CRISPOR uses BWA, a few tools from the UCSC Genome Browser (twoBitToFa, bedClip),
various R packages and a huge collection of external packages and source code files
from published articles, see the file crisporEffScores.py for the exact references.

Installation of the package:

    make crispor_env

    source crispor_env/bin/activate

    make devbuild

Install required R libraries:
   
    sudo Rscript -e 'install.packages(c("e1071"),  repos="http://cran.rstudio.com/")'
    sudo Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("limma"));'

Now in a python console type this:

    from crispor_cli import crispor
    
    crispor.main(args,options)

Description for args and options is given below:
  args=[<genome_name>,<input_fasta>,<output_file>]
Here args is a list containing org,fastaInFile and guideOutFile in this order:
  Example args -
    args=['sacCer3','/input/guide_yeast.fasta','/output/yo_guide.tsv']

And options is a dictionary containing all the extra options permitted by crispor.
  Example options - 
    options = {'offtargetFname':'/output/yo_off.tsv','pam':'NGG','debug':True,'skipAlign':True}

Here are the keys that can be added to options dictionary-

Options:
  debug      -     show debug messages, do not delete temp directory
  test      -      run internal tests
  pam    -              PAM-motif to use, default NGG. TTTN triggers special
                        Cpf1 behavior: no scores anymore + the PAM is assumed
                        to be 5' of the guide. Common PAMs are:
                        NGG,TTTN,NGA,NGCG,NNAGAA,NGGNG,NNGRRT,NNNNGMTT,NNNNACA
  offtargetFname - 
                        write offtarget info to this filename
  maxOcc - 
                        MAXOCC parameter, guides with more matches are
                        excluded

  mismatches-
                         maximum number of mismatches, default 4
  
  skipAlign  -
                        do not align the input sequence. The on-target will be
                        a random match with 0 mismatches.
  noEffScores -
                        do not calculate the efficiency scores
  minAltPamScore -
                        minimum MIT off-target score for alternative PAMs, default
                        1.0
  genomeDir-
                        directory with genomes, default ./genomes
```
    

# Licenses

Included software:

* BWA is under GPL3
* libSVM: under copyright by Chih-Chung Chang and Chih-Jen Lin see http://www.csie.ntu.edu.tw/~cjlin/libsvm/COPYRIGHT
* svmlight: free for non-commercial use, see http://svmlight.joachims.org/
* SSC: no license specified
* primer3: GPL2.
* Fusi/Doench score: see LICENSE.txt, (c) by Microsoft Research
* crispor.py and crisporEffScores.py themselves are released under GPLv3, see LICENSE.txt
