
```
GitHub Conventions:
   - Master branch is the current release
   - Development branch is not gauranteed to be stable
   - Source releases are avaialable under the GitHub release tab and at 
     http://www.repeatmasker.org
```

RepeatModeler
=============

RepeatModeler is a de novo transposable element (TE) family identification and 
modeling package. At the heart of RepeatModeler are three de-novo repeat finding
programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ 
complementary computational methods for identifying repeat element boundaries
and family relationships from sequence data.

RepeatModeler assists in automating the runs of the various algorithms
given a genomic database, clustering redundant results, refining and 
classifying the families and producing a high quality library of
TE families suitable for use with RepeatMasker and ultimately for submission
to the Dfam database ( http://dfam.org ).

Authors
-------
 RepeatModeler:
   Robert Hubley, Arian Smit - Institute for Systems Biology

 LTR Pipeline Extensions:
   Jullien M. Flynn - Cornell University

Installation
------------

 There are two supported paths to installing RepeatModeler on a 
 UNIX-based server. RepeatModeler may be installed from source as
 described in the "Source Distribution Installation" instructions
 below, or using one of our Dfam-TETools container images ( Docker or 
 Singularity ). The containers include RepeatModeler, it's 
 prerequisites and additional TE analysis tools/utilities used by
 Dfam. Information on the Dfam-TETools container may be found here:
 https://github.com/Dfam-consortium/TETools
 

Source Distribution Installation
--------------------------------

**Prerequisites**

  Perl
    Available at http://www.perl.org/get.html. Developed and tested
    with version 5.8.8.

  RepeatMasker & Libraries
    Developed and tested with open-4.1.4. The program is available at 
    http://www.repeatmasker.org/RMDownload.html and is distributed with
    Dfam - an open database of transposable element families.

  RECON - De Novo Repeat Finder, Bao Z. and Eddy S.R.
    Developed and tested with our patched version of RECON ( 1.08 ). 
    The 1.08 version fixes problems with running RECON on 64 bit machines and 
    supplies a workaround to a division by zero bug along with some buffer 
    overrun fixes. The program is available at: 
       http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz. 
    The original version is available at http://eddylab.org/software/recon/.

  RepeatScout - De Novo Repeat Finder, Price A.L., Jones N.C. and Pevzner P.A.
    Developed and tested with our multiple sequence version of RepeatScout 
    ( 1.0.6 ).  This version is available at 
    http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz

  TRF - Tandem Repeat Finder, G. Benson et al.
    You can obtain a free copy at http://tandem.bu.edu/trf/trf.html
    RepeatModeler requires version 4.0.9 or higher.

  RMBlast - A modified version of NCBI Blast for use with RepeatMasker
    and RepeatModeler.  Precompiled binaries and source can be found at
    http://www.repeatmasker.org/RMBlast.html
    We highly recommend using 2.13.0 or higher.

  UCSC genome browser command-line utilities - Some tools included with
    RepeatModeler work with files in the 'twobit' file format using
    these programs: twoBitToFa, faToTwoBit, and twoBitInfo.
    Precompiled binaries and source for these programs can be found at
    http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads.

**Optional. Required for running LTR structural search pipeline:**

  LtrHarvest - The LtrHarvest program is part of the GenomeTools suite.  We
    have developed this release of RepeatModeler on GenomeTools version 1.5.9 
    available for download from here: http://genometools.org/pub/
    NOTE: use the "make threads=yes" build options to enable multi-threaded
          runs.
          
  Ltr_retriever - A LTR discovery post-processing and filtering tool.  We 
    recommend using version 2.6 or higher from here: 
    https://github.com/oushujun/LTR_retriever/releases

  MAFFT - A multiple sequence alignment program.  We developed and tested
    RepeatModeler using mafft version 7.407.  Please use this verison or
    higher from here:
    https://mafft.cbrc.jp/alignment/software/

  CD-HIT - A sequence clustering package.  We developed and tested
    RepeatModeler using version 4.8.1.  Please use this version or higher
    from:
    http://weizhongli-lab.org/cd-hit/

  Ninja - A tool for large-scale neighbor-joining phylogeny inference and 
    clustering.  We developed and tested RepeatModeler using Ninja version
    "0.98-cluster_only". Please obtain a copy from:
    https://github.com/TravisWheelerLab/NINJA/releases/tag/0.98-cluster_only

**Installation**

  1. Obtain the source distribution
       - Github : https://github.com/Dfam-consortium/RepeatModeler
         Available by cloning the master branch of the RepeatModeler repository
         ( latest released version ) or by downloading a release from the 
         repository "releases" tab.
     or 
       - RepeatMasker Website : http://www.repeatmasker.org/RepeatModeler

  2. Uncompress and expand the distribution archive:

     ```
     tar -zxvf RepeatModeler-#.#.tar.gz
     ```

  3. Configure for your site:

     Automatic:
       + Run the "configure" script interactively with prompts
         for each setting:

               perl ./configure

       + Run the "configure" script with supplied paramters:

               perl ./configure -rscout_dir .. -recon_dir ..
 
     By Hand:
       + Edit the configuration file "RepModelConfig.pm"

     Dynamically:
       + Use the "configuration overrides" command line options
          with the RepeatModeler programs. e.g:

               ./RepeatModeler -rscout_dir .. -recon_dir .. 
         
 
Example Run
-----------

In this example we first downloaded elephant sequences
from Genbank ( approx 11MB ) into a file called elephant.fa. 

  1. Create a Database for RepeatModeler

     RepeatModeler uses a NCBI BLASTDB as input to the
     repeat modeling pipeline.  A utility is provided to assist
     the user in creating a single database from several 
     types of input structures.  

           <RepeatModelerPath>/BuildDatabase -name elephant elephant.fa

     Run "BuildDatabase" without any options in order to see the 
     full documentation on this utility. There are several options
     which make it easier to import multiple sequence files into
     one database.
  
     TIP: It is a good idea to place your datafiles and run this
          program suite from a local disk rather than over NFS. 
          This will greatly improve runtime as the filesystem 
          access is considerable

  2. Run RepeatModeler
     
     RepeatModeler runs several compute intensive programs on the
     input sequence.  For best results run this on a single machine with
     a moderate amount of memory > 32GB and multiple processors.  
     Our setup is Xeon(R) CPU E5-2680 v4 @ 2.40GHz - 28 cores, 128GB RAM.
     To specify a run using 20 threads (at most), and including the new 
     LTR discovery pipeline:
     
            nohup <RepeatModelerPath>/RepeatModeler -database elephant 
            -threads 20 -LTRStruct >& run.out &

     The nohup (or screen) is used on our machines when running long 
     jobs.  The log output is saved to a file and the process is backgrounded.
     For typical runtimes ( can be > 1-2 days with this configuration on a
     well assembled mammalian genome ) see the run statistics section of 
     this file.  It is important to save the log output for later usage.  
     It contains the random number generator seed so that the sampling 
     process may be reproduced if necessary.  In addition the log file 
     contains details about the progress of the run for later assesment 
     of peformance or debuging problems.

  3. Interpret the results

     RepeatModeler produces a voluminous amount of temporary files stored
     in a directory created at runtime named like:

               RM_<PID>.<DATE> ie. "RM_5098.MonMar141305172005" 

     and remains after each run for debugging purposes or for the purpose
     of resuming runs if a failure occures. At the succesful completion 
     of a run, three files are generated:

               <database_name>-families.fa  : Consensus sequences
               <database_name>-families.stk : Seed alignments
               <database_name>-rmod.log     : A summarized log of the run

     The seed alignment file is in a Dfam compatible Stockholm format and
     may be uploaded to the Dfam database by submiting the data to 
     help@dfam.org or by going to dfam.org/login and creating an upload
     account.

     The fasta format is useful for running quick custom library searches
     using RepeatMasker.  Ie.:

               <RepeatMaskerPath>/RepeatMasker -lib <database_name>-families.fa mySequence.fa 

     Other files produced in the working directory include:

```
       RM_<PID>.<DATE>/
          consensi.fa
          families.stk
          round-1/
               sampleDB-#.fa       : The genomic sample used in this round
               sampleDB-#.fa.lfreq : The RepeatScout lmer table
               sampleDB-#.fa.rscons: The RepeatScout generated consensi
               sampleDB-#.fa.rscons.filtered : The simple repeat/low 
                                               complexity filtered
                                               version of *.rscons
               consensi.fa         : The final consensi db for this round
               family-#-cons.html  : A visualization of the model
                                     refinement process.  This can be opened
                                     in web browsers that support zooming.
                                     ( such as firefox ).
                                     This is used to track down problems 
                                     with the Refiner.pl
               index.html          : A HTML index to all the family-#-cons.html
                                     files.
          round-2/
               sampleDB-#.fa       : The genomic sample used in this round
               msps.out            : The output of the sample all-vs-all 
                                     comparison
               summary/            : The RECON output directory
                    eles           : The RECON family output
               consensi.fa         : Same as above
               family-#-cons.html  : Same as above
               index.html          : Same as above
          round-3/
               Same as round-2
           ..
          round-n/
```


  Recover from a failure

  If for some reason RepeatModeler fails, you may restart an
  analysis starting from the last round it was working on.  The
  -recoverDir [<i>ResultDir</i>] option allows you to specify a
  diretory ( i.e RM_<PID>.<DATE>/ ) where a previous run of
  RepeatModeler was working and it will automatically determine
  how to continue the analysis.

Caveats
-------

  + RMBlast uses the NCBI stat reporting mechanism to report
    usage statistics back over the net.  If RepeatModeler is 
    taking in inordinate amount of time to complete 
    ( > 1 week on a multi-core machine ) or you do not have
    an outside network connection to the machine, you should
    disable this reporting feature by setting the environment
    variable BLAST_USAGE_REPORT=false or by creating a .ncbirc
    file in the users home directory with the stanza:

    [BLAST]

    BLAST_USAGE_REPORT=false

  + RepeatModeler is designed to run on assemblies rather
    than genome reads. At the start of a run a quick analysis
    is performed on the input database to ascertain the 
    assembly N50.  A histogram of contig size is also displayed.

  + RepeatModeler employs symmetric multiprocessing parallelism,
    therefore should be run on a single machine per-assembly.

  + It is not recommended that a genome be run in a batched fashion
    nor the results of multiple RepeatModeler runs on the same 
    genome be naively combined.  Doing so will generate a combined
    library that is largely redundant.  The -genomeSampleSizeMax 
    parameter is provided for the purpose of increasing the amount
    of the genome sampled while avoiding rediscovery of families.

  Please see the RELEASE-NOTES file for more details.


RepeatModeler Statistics
------------------------
Benchmarks and statistics for runs of RepeatModeler on several sample
genomes.

RepeatModeler 2.0 ( RECON + RepeatScout + LTRStruct ):

Genome          |  Genome DB Size (bp) |   Runtime (hh:mm)*  |  Models Built
----------------|----------------------|---------------------|--------------
D. melanogaster |  164 Mbp             |    12:56            |     734   
D. rerio        |  1.4 Gbp             |    40:36            |    3851   
O. sativa       |  375 Mbp             |    37:23            |    2648   

  *  Analysis run on a CentOS 7.6.1810 Linux system with
     Intel... processors.  D. melanogaster was run using 
     8 parallel jobs ( -pa 8 ) while D. rerio and O.sativa
     were run with 16 parallel jobs.

Previous Versions

```
RepeatModeler 1.0.3 ( RECON + RepeatScout ):

            Genome DB    Sample***   Run Time*   Models   Models       % Sample
Genome      Size (bp)    Size (bp)   (hh:mm)     Built    Classified   Masked**
----------  -----------  ----------  ----------  -------  -----------  --------
Marmoset    3.0 Bbp      228 Mbp      55:37        582      564         36.4    
Zebrafinch  1.2 Bbp      222 Mbp      66:29        178       75          8.6    


  *  Analysis run on a 4 processor P4, 2.4Ghz, 3GB RAM, machine
     running Red Hat Linux.

  ** Includes simple repeats and low complexity DNA. Results
     obtained with RepeatMasker open-3.2.5, WUBlast and
     the -lib option.

  *** Sample size does not include 40 Mbp used in the RepeatScout analysis.  
      This 40 Mbp is randomly chosen and may overlap 0-100% of the 
      sample used in the RECON analysis.


RepeatModeler 1.0.2 ( RECON + RepeatScout ):

            Genome DB    Sample***   Run Time*   Models   Models       % Sample
Genome      Size (bp)    Size (bp)   (hh:mm)     Built    Classified   Masked**
----------  -----------  ----------  ----------  -------  -----------  --------
Human HG18     3.1 Bbp      238 Mbp   46:36        614      611         35.66
Platypus       1.9 Bbp      220 Mbp   76:02        540      457         
Zebrafinch     1.3 Bbp      220 Mbp   63:57        233      104          9.41
Sea Urchin     867 Mbp      220 Mbp   40:03       1830      360         33.85
diatom       32,930,227  32,930,227    4:41        128       35          2.86
Rabbit       11,770,949  11,770,949    3:14         83       72         31.30

  *  Analysis run on a 4 processor P4, 2.4Ghz, 3GB RAM, machine
     running Red Hat Linux.

  ** Includes simple repeats and low complexity DNA. Results
     obtained with RepeatMasker open-3.1.9, WUBlast and
     the -lib option.

  *** Sample size does not include 40 Mbp used in the RepeatScout analysis.  
      This 40 Mbp is randomly chosen and may overlap 0-100% of the 
      sample used in the RECON analysis.



RepeatModeler 1.0.0 ( RECON ):

            Genome DB    Sample      Run Time*   Models   Models       % Sample
Genome      Size (bp)    Size (bp)   (hh:mm)     Built    Classified   Masked**
----------  -----------  ----------  ----------  -------  -----------  --------
Opossum     3.5 Billion  319 Mbp     140:52      1137     467          52.55
Armadillo   47,332,136   47,332,136    6:20      121      92           36.07
Platypus    14,768,992   14,768,992    3:46      18       13           40.69
Rabbit      11,770,949   11,770,949    2:17      20       16           28.67
Elephant    11,550,090   11,550,090    1:21      34       28           37.08

  *  Analysis run on a 4 processor P4, 2.4Ghz, 3GB RAM, machine
     running Red Hat Linux.
  ** Includes simple repeats and low complexity DNA. Results
     obtained with RepeatMasker open-3.0.9, WUBlast and
     the -lib option.
```



Credits
-------

Arnie Kas for the work done on the original MultAln.pm.

Andy Siegel for statistics consultations.

Thanks so much to Warren Gish for his invaluable assistance 
and consultation on his ABBlast program suite.

Alkes Price and Pavel Pevzner for assistance with RepeatScout
and hosting my multi-sequence version of RepeatScout.

Shujun Ou, and Ning Jiang for discussions and assistance with 
using LTR_retreiver.

This work was supported by the NIH ( R44 HG02244-02), 
( RO1 HG002939 ), ( U24 HG010136 ), and the Institute 
for Systems Biology.


License
-------

This work is licensed under the Open Source License v2.1.  
To view a copy of this license, visit 
http://www.opensource.org/licenses/osl-2.1.php or
see the LICENSE file contained in this distribution.

