#Running TULIP (The Uncorrected Long-red Integration Process), version 0.4 late 2016 (European eel)

TULIP currently consists of to Perl scripts, tulipseed.perl and tulipbulb.perl. These are very much intended as prototypes, and additional components and/or implementations are likely to follow. <br>
Tulipseed takes as input alignments files of long reads to sparse short seeds, and outputs a graph and scaffold structures. Tulipbulb adds long read sequencing data to these.<br>

##The steps below describe how we used TULIP to assemble the European eel genome:<br>


1. **Input data**


  _eel_seeds_285.fasta:_
  These are pre-selected seed sequenced, with criteria 'not too repetitive'. Seeds containing repetitive sequence should be fine, but it will then take much longer to untangle the graph, and to optimize seed numbers. Strict requirements are that seed sequence names are unique numbers (excluding 0), and that all seeds are of exactly the same length.<br>
  ```
Example file:
>53
TCTGTATAATTTTTTTTTTTTCAAAAAGAAATGTCACAATTATTCACACCCCAGTTTTCAGCACCCTCTCTTAACAAGGACATTCTTCTGTAATATTTTATGAGATAGATGGACACATCCTTGTCCATTCTTGCATACACCATCTTTCTAAATTTTCTACTGAAAATGTCCTCCTCAGTTCAAACCAGAAAATTTGGTTACATTCTGGAAACTTGAATATTGATCCAGAGACAAAAACAGCAAAACAGTAATTTTGTGGTAAATTAATCATTTATTGGTTGAGTT  
>118
GTTGCAAGCATATTTTAGCATTCCTTTAGCTCAAAAGTTTCTCATTTTTTTCTTGCCCATTATCAACAGTGACAAATTCTTCTGATATACATCTTTCTGATGTTTGTGGTTCCACATTGGCCTTCTCCTGCATTGTGGTATTTCTACTTTGTTTAGTTAATCAGCTGTTGAAATTAGCCTTTAGTCCCACAGGGAATTACAGGAATTGTGGTATACACTGTTATAAGCAATATACATTTTATTTTATGATACCTGCTAAAGAAGGTAATATGTCAGATGTTATAG  
>154
TTCTGAATTCCTTTAAGACTTCAAGGTGAATGGTGAATTAAAGTGCTGCCATCATATAGGCTGTTTAAAGGCAGTTTTAAATGATTTTATATATATTTTATATGATTACAGACAATGTGATTCATGAAGAAAATGTGGGCAGTCCTTTTCCCTGTAGCAAGGTCAGTAAAATAATAGTGACAGAATAATGTGCTTGACGTCTCTAATTTTACAATCTCATATACCACTGTATGCCTATGTGAGTCAAATATGATATAAAATTGAACATTATTATGTTTGTAATGG  
  ```
  _long_reads.fasta:_
  Long reads in FASTA format. The read names should be unique in each file, so be careful if you concatenate data!<br><br>


2. **Alignments**


  We used BWA MEM to align reads to seeds, but other aligners should work. Currently, TULIP accepts SAM format (support for other formats, e.g. DALIGNER output, was present in even earlier versions and might reappear). Only the first 6 fileds of the SAM alignment information is used (up to and including the CIGAR string), so you might want to clip off the rest.

  Command lines:
  ```
  bwa index -p 285_seeds eel_seeds_285.fasta
  bwa mem -t 4 -k 14 -W 45 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 285_seeds R7.3.fasta			| cut -f 1,2,3,4,5,6 > R73_vs_285.shortsam  
  bwa mem -t 4 -k 16 -W 50 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 285_seeds R9_pass_1d.fasta	| cut -f 1,2,3,4,5,6 > R91D_vs_285.shortsam  
  bwa mem -t 4 -k 19 -W 60 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 285_seeds R9_pass_2d.fasta	| cut -f 1,2,3,4,5,6 > R92D_vs_285.shortsam  
  bwa mem -t 4 -k 16 -W 60 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 285_seeds R9.4_pass_1d.fasta	| cut -f 1,2,3,4,5,6 > R941D_vs_285.shortsam  
  ```
<br>


3. **Configuration**


  When using multiple alignment files, you should specify their locations in a short configuration file, format type [tab] alignment [tab] original fasta:
  ```
  sam	R92D_vs_285.shortsam	R9_pass_2d.fasta
  sam	R91D_vs_285.shortsam	R9_pass_1d.fasta
  sam	R941D_vs_285.shortsam	R9.4_pass_1d.fasta
  sam	R73_vs_285.shortsam		R7.3.fasta
  ```
<br>

4. **TULIP seed layout**

  ```
  ./tulipseed.perl --seedlength 285 --config alignments.txt --diploid --out tulip/eel
  ```
  
  _This will generate the following files:_<br>
  ```
  tulip/eel.graph				The seed graph, text format
  tulip/eel.graph_tmp			The seed graph, binary format
  tulip/eel.scaffolds			Ordered seed scaffolds, text format
  tulip/eel.scaffolds_tmp		Ordered seed scaffolds, binary format
  tulip/eel.scaffolds			Ordered seed scaffolds, text format
  tulip/eel.seeds_tmp			Seed usage information
  tulip/eel.layout_log			A log file
  tulip/eel.scaffolds_stats		Length statistics per scaffold
  ```
  The _.*_tmp_ files will be used by tulipbulb.perl.<br>
  The _.graph_ file simply list which seeds are connected in the final simplified graph:<br>
  
  ```
  100002	ii	3501405	10		1		 432	 559.50		 585	 43.49	2284
  100002	oo	5662469	12		1		 -60	 -56.92		 -53	  2.75	2284
  1000042	ii	2229111	 5		1		2120	2370.40		2587	190.04	1019
  1000042	oi	6256249	 4		1		3531	3656.50		3963	178.04	1019
  1000048	io	5919951	 6		1		1509	1540.33		1628	 40.46	1598
  1000048	oi	672698	 7		1		1107	1149.57		1225	 34.42	1598
  1000049	ii	239093	 3		1		 629	 646.33		 668	 16.21	1193
  ```
  
  
  _Columns indicate:_<br>
  ```
  * Seed 1			Name of the first seed
  * Orientation		ii | io | oi | oo, how seeds are connected (in/out)
  * Seed 2			Name of the second seed
  * Evidence		The number of long read alignments connection these seeds
  * Hypothetical	1 = actual link, 0 = inferred link
  * Minimum			Minimum gap between the seeds observed in the alignments
  * Mean			Mean gap between the seeds observed in the alignments
  * Maximum			Maximum gap between the seeds observed in the alignments
  * StDev			Standard deviation of the gap estimate
  * Scaffold		Final scaffold ID
  ```
<br>  

5. **TULIP bundling**
  ```
  ./tulipbulb --seeds eel_seeds_285.fasta --config alignments.txt --input tulip/eel --out bulb/eel
  ```

  This will add sequence from the original seeds and reads to the graph output by tulipseed.perl.

  Output files are:
  ```
  bulb/eel_readbundle_999.fasta		The reads used to construct scaffold 999
  bulb/eel_scaffold_999.fasta		The sequence for scaffolds 999
  bulb/eel.bundle_log				A log file
  ```
  
  For each scaffold, two files are generated. The scaffold sequence shows sequence derived from seeds and long reads in upper and lower case, respectively.<br><br>

##N.B. Prototype could include N's
The tulipbulb.perl script contains a bug which prevents it from adding sequence in rare cases (it will then add gaps, NNNs). We will fix this soon, and also add some additional output options.<br>

##For any questions please send us an e-mail using:
```
    m.liem@biology.leidenuniv.nl
c.v.henkel@biology.leidenuniv.nl
```
