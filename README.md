#HUPAN

---

 **1. Introduction**
The human reference genome is still incomplete, especially for those population-specific or individual-specific regions, which may have important functions. It encourages us to build the pan-genome of human population. Previously, our team developed a "map-to-pan" strategy--[EUPAN][1], specific for eukaryotic pan-genome analysis. However, due to the large genome size of individual human genome, [EUPAN][2] is not suit for pan-genome analysis involving in hundreds of individual genomes. Here, we present an improved tool, HUPAN (Human Pan-genome Analysis), for human pan-genome analysis.

**2. Installation**

**Requirements** 

 - R 3.1 or later (https://www.r-project.org/)
    
    R is utilized for visualization and statistical tests in HUPAN
    toolbox. Please install R first and make sure R and Rscript are
    under your PATH.
 - R packages Several R packages are needed including ggplot2, reshape2
    and ape packages. Follow the Installation step,
 - or you can install the packages by yourself.

**Installation procedures** 

 - Download the HUPAN toolbox from [github][3]:

    git clone git@github.com:SJTU-CGM/HUPAN.git

 - Alternatively, you also could obtain the package in the
   http://cgm.sjtu.edu.cn/hupan/. Please uncompressed the HUPAN toolbox
   package.

    tar zxvf HUPAN-v**.tar.gz 

 - Install necessary R packages

    cd HUPAN & Rscript installRPac 

 - Compile necessary tools.
   
   make
   
     You will find executable files: *ccov*, *bam2cov* and *hupan* et
   al. in bin/ directory.

 - Add bin/ to PATH and add lib/ to LD_LIBRARY_PATH. To do this, add the
   following text to ~/.bash_profile
   

       export PATH=$PATH:/path/to/HUPAN/bin:   
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/HUPAN/lib/:   
       export PERL5LIB=$PERL5LIB:/path/to/HUPAN/lib/:

 

 - and run 
   
       source ~/.bash_profile


 - Test if HUPAN toolbox is installed successfully. `hupan` If you see
   the following content, congratulations! HUPAN toolbox is successfully
   installed. If not, see if all the requirements are satisfied or
   contact the author for help.

   

     Usage: hupan <command> ...
        
        Available commands:
        	qualSta      	View the overall sequencing quality of a large number of files
        	trim         	Trim or filter low-quality reads parallelly
        	alignRead    	Map reads to a reference parallelly
        	sam2bam      	Covert alignments (.sam) to sorted .bam files
        	bamSta       	Statistics of parallel mapping
        	assemble     	Assemble reads parallelly
        	alignContig  	Align assembly results to a referenece parallelly
        	extractSeq   	Extract contigs parallelly
        	assemSta     	Statistics of parallel assembly
        	getUnalnCtg  	Extract the unaligned contigs from nucmer alignment (processed by quast)
        	rmRedundant  	Remove redundant contigs of a fasta file
        	pTpG         	Get the longest transcripts to represent genes
        	geneCov      	Calculate gene body coverage and CDS coverage
        	geneExist    	Determine gene presence-absence based on gene body coverage and CDS coverage
        	subSample    	Select subset of samples from gene PAV profile
        	gFamExist    	Determine gene family presence-absence based on gene presence-absence
        	bam2bed      	Calculate genome region presence-absence from .bam
        	fastaSta     	Calculate statistics of fasta file
        	sim          	Simulation and plot of the pan-genome and the core genome
        	getTaxClass  	Obtain the taxonomic classification of sequences
        	rmCtm        	Detect and discard the potentail contamination
        	blastAlign   	Align sequences to target sequence by blast
        	simSeq       	Simulate and plot the total size of novel sequences
        	splitSeq     	Split sequence file into multiple small size files
        	genePre      	Ab initio gene predict combining with RNA and protein evidence
        	mergeNovGene 	Merge maker result from multiple maker result files
        	filterNovGene	Filter the novel precited genes.

**3.	Main analysis procedures**
HUPAN is an improved version of EUPAN, and most functions were directly obtained from EUPAN. Besides, several new functions were developed for efficient pan-genome analysis on human genome. Here, we provide the main analysis procedures of human pan-genome analysis on an example data. And we provide three types of tools:

 1. Single machine version;
 2. LSF version (working on supercomputer based on LSF system, in
    which, "bsub" is used to submit jobs);
 3. SLURM version (working on supercomputer based on SLURM system, in
    which, "sbatch" is used to submit jobs).

Due to the large genome size of individual genome, conducting pan-genome analysis on hudreds of individuals could hardly finish on in single machine. We strongly suggestted the users conduct all the analysis in the supercomputer implemented LSF system or SLURM system. All the commands of `hupanLSF` and `hupanSLURM` are same excepted for the way of submit jobs are different. In the following, we give all the exmaples of command based on SLURM system. If the users work on supercomputer based on LSF system, please replace “`hupanSLURM`” with “`hupanLSF`”.
**(1)	Example data**
This data set includes sequencing data of three samples from [NA12878][4]. Each sample include 6,000,000 paired-end reads that could map to chromosome 22. Note these are only simple example to help users understand the input data type and data structure, and run the pipeline. The real data may be much larger and more complex. 
Please download [here][5] and undecompress it:

    tar zxvf hupanExample.tar.gz & cd hupanExample

And you can find two directories:

    data/   The sequencing data with a per-sample-per-directory structure;
    ref/    The reference sequence and annotation information (chr22 of GRCh38).

**(2) The parallel quality control**
The tools [FastQC][6] and [Trimmomatic][7] were used to view the sequencing quality and trim low-quality reads.
i. The command `qualitySta` is used to overview qualities:

    hupanSLURM  qualSta -f /path/to/Fastqc -t 16 -v PE data/ preview_quality/

Results can be found in the `preview_quality/` directory.
ii. If the reads are not so good, the users could trim or filter low-quality reads by the command `trim`:

    hupanSLURM trim data/ trim/ /path/to/Trimmomatic
    hupanSLURM trim -w 100 -m 100 data/ filter/ /path/to/Trimmomatic

Results could be found in the trim or filter directory.
iii.After trimming or filtration of reads, the sequencing quality should be evaluated again by `qualitySta`, and if the trimming results are still not good for subsequent analyses, new parameters should be given and the above steps should be conducted for several times.
**(3) De novo assembly of individual genomes**
To obtain non-reference sequences from each individual genome, we need first to conduct de novo assembly on the raw reads. We provide three distinct strategies:
i.Directly assembly by [SOAPDenovo2][8]. Please note that this startegy requires huge memory for assembly an individual human genome (according to our test, finishing the assembly of a human genome of 30-fold sequencing data needs more than 500 Gb memory), we strongly suggested that do not use this command unless you have multiple supercomputer with huge memory.

    hupanSLURM assemble soapdenovo -t 16 -k 91 data/ assembly_soap/ /path/to/SOAPDenovo2/

ii.Assembly by the iterative use of SOAPDenovo2. Not Recommended.

    hupanSLURM linearK data assembly_linearK/ /path/to/SOAPDenovo2 

iii. Assembly by [SGA][9]. We recommend the users preform assembly by this command. According to our experience on 185 newly sequenced genomes, the maximum memory consumption in assembling the human genome of 30-fold sequencing data is about 60Gb.

    hupanSLURM assemble sga -t 16 data/ assembly_result /path/to/sga/

**(4) Extract non-reference sequences from assembled contigs** 	
i. In order to obtain the non-reference sequence from each individual genome, the assembled contigs are aligned to the reference genome with nucmer tool within [Mummer][10] package.

    hupanSLURM alignContig assembly_result/data/ aligned_result	 /path/to/MUMmer/ /path/to/reference.fa

ii. Then the contigs those are highly similar with the reference genome are discarded and the remaining contigs are considered as candidate non-reference sequences.

    hupanSLURM extractSeq assembly_result/data/ candidate aligned_result

iii. All the candidate non-reference sequences are assessed by [QUAST][11] to obtain non-reference sequences.  

    hupanSLURM assemSta candidate/data/ quast_result /path/to/quast-4.5/ /path/to/reference.fa

iv. Two types of non-reference sequences, fully unaligned sequences and partially unaligned sequences for each individual could be collected:

    hupanSLURM getUnalnCtg -p .contig candidate/data/ quast_result/data/ Unalign_result

v. Non-reference sequences from multiple individuals are merged.

    hupanSLURM mergeUnalnCtg Unalign_result/data/ mergeUnalnCtg_result

**(5) Remove redundancy and potential commination sequences**
i. After obtaining the non-reference sequences from multiple individuals, redundant sequences between different individuals should be excluded, and the potential commination sequences from non-human species are also removed for further analysis. The step of remove redundancy sequences is conducted by [CDHIT][12] for fully unaligned sequences and partially unaligned sequences, respectively.

    hupanSLURM rmRedundant cdhitCluster  mergeUnalnCtg_result/total.fully.fa rmRedundant.fully.unaligned /path/to/cdhit/
    hupanSLURM rmRedundant cdhitCluster mergeUnalnCtg_result/total.partilly.fa rmRedundant.partially.unaligned /path/to/cd-hit/

ii. Then the non-redundant sequences are aligned to NCBI’s non-redundant nucleotide database by blastn. 

    hupanSLURM blastAlign blast rmRedundant rmRedundant_blast /path/to/nt /path/to/blast

iii. According to the alignment result, the taxonomic classification of each sequences (if have) could be obtained.

    hupanSLURM getTaxClass rmRedundant_blast/ data/fully/fully.non-redundant.blast info/ TaxClass_fully
    hupanSLURM getTaxClass rmRedundant_blast/ data/partially/partially.non-redundant.blast info/ TaxClass_partially

iv. And the sequences classifying as microbiology and non-primate eukaryotes are considered as non-human sequences and removed from further consideration  

    hupanSLURM rmCtm -i 60 rmRedundant/fully/fully.non-redundant.fa rmRedundant_blast/data/fully/fully.non-redundant.blast TaxClass_fully/data/accession.name rmCtm_fully
    hupanSLURM rmCtm -i 60 rmRedundant/partially/partially.non-redundant.fa rmRedundant_blast/data/partially/partially.non-redundant.blast TaxClass_partially/data/accession.name rmCtm_partially

**(6) Construction and annotation of pan-genome**
i. The non-redundant sequences of fully unaligned sequences and partially unaligned sequences are merged and further clustered to remove redundant sequences.

    mkdir Nonreference
    cat rmCtm_fully/data/novel_sequence.fa rmCtm_partially/data/novel_sequence.fa > Nonreference/nonrefernce.before.fa
    hupanSLURM rmRedundant cdhitCluster Nonreference/nonrefernce.before.fa NonredundantNonreference /path/to/cdhit/

ii. And the resulted sequences together with the human reference genome construct the pan-genome sequences. The annotation of reference genome could be directly download from [GenCODE][13] or other common used annotation dataset. The annotation information of non-reference sequences is predicted by [MAKER][14].Usually, the size of non-reference sequences is large and the procedure of gene prediction is slow. We recommend the users to split the file of non-reference sequences into multiple small files and run maker in parallel.

    hupanSLURM splitSeq NonredundantNonreference/non-redundant.fa GenePre_input
    hupanSLURM genePre GenePre_input GenePre /path/to/maker/config_file /path/to/maker

iii. Then after all procedures are finished, the outcomes are merged.

    hupanSLURM mergeNovGene GenePre GenePre_merge /path/to/maker

iv. The new predicted gene may be highly similar with the reference genome, and additional filtering step should be conducted to ensure the novelty of predicted gene. 

    hupanSLURM filterNovGen GenePre_merge GenePre_filter /path/to/reference/ /path/to/blast /path/to/cdhit /path/to/RepeatMask

v. The annotation of pan-genome sequences is simply to obtain by combine two annotation file.
 

     hupanSLURM pTpG ref/ref.gtf ref/ref-ptpg.gtf
      cat ref/ref-ptpg.gtf non-reference.gtf >pan/pan.gtf

**(7) PAV analysis**
The “map-to_pan” strategy is utilized to determine gene presence-absence. 
i. The raw reads are mapped to pan-genome sequences by [Bowtie2][15]

      cd pan & /path/to/bowtie2/bowtie2-build pan.fa pan &cd ..
      hupanSLURM alignRead –f bowtie2 data/ map2pam /path/to/bowtie2 pan/pan

ii. The result of .sam should be converted to .bam and sorted and indexed use [Samtools][16].
  

    hupanSLURM sam2bam map2pan/data panBam /path/to/samtools

iii. Then the gene body coverage and the cds coverage of each gene are calculated.

      hupanSLURM geneCov panBam/data geneCov/ pan/pan.gtf

iv. Finally, the gene presence-absence is determined by the threshold of cds coverage as 95%.

      mkdir geneExist & hupanSLURM geneExist geneCov/summary_gene.cov geneCov/summary_cds.cov 0 0.95 >geneExist/gene.exist

**(8) Bugs or suggestions**
Any bugs or suggestions, please contact the [authors][17]. 



 


  [1]: http://cgm.sjtu.edu.cn/eupan/
  [2]: http://cgm.sjtu.edu.cn/eupan/
  [3]: https://github.com/SJTU-CGM/HUPAN
  [4]: ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
  [5]: http://cgm.sjtu.edu.cn/hupan/data/ExampleData.tar.gz
  [6]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  [7]: http://www.usadellab.org/cms/index.php?page=trimmomatic
  [8]: https://sourceforge.net/projects/soapdenovo2/
  [9]: https://github.com/jts/sga/
  [10]: https://sourceforge.net/projects/mummer/
  [11]: http://quast.bioinf.spbau.ru/
  [12]: http://weizhongli-lab.org/cd-hit/
  [13]: https://www.gencodegenes.org/
  [14]: http://www.yandell-lab.org/software/maker.html
  [15]: http://bowtie-bio.sourceforge.net
  [16]: http://www.htslib.org/
  [17]: zhqduan@sjtu.edu.cn