# ChipseqNewdataPipeline_2
A bash pipeline for ChIP - seq data processing

Welcome to ChipseqNewdataPipeline_2!

We hope this program will help you process and analyze your ChIP-seq data. It contains a series of scripts that will first create peak expression files and will then run an R script to give information about your transcriptor factor's target genes and their functions. Finally, the HOMER findMotifsGenome tool will be used to define the sequence that binds DNA.

To make sure the program works properly, prior to launching it, you need to make sure your genome file is unzipped, and install the following:

1. The HOMER program and the information of the organism of study (in the test file of *Arabidopsis* transcriptor factor it's tair10)
	- Please, make sure you add your script file to your PATH so that HOMER can produce the desired results.
2. The following R packages:
	- ChIPseeker package 
	- TxDB of the organism of study
	- .db file of the organism of study (org.At.tair.db for *Arabidopsis*)
	- clusterProfiler
	- pathview

To run the program, please make sure you filla ll the test_params.txt parameters correctly and then launch both with the ./ command. Whenever you launch ChipseqNewdataPipeline without parameters, it will also ask for this information.

If you have any further questions, do not hesitate to ask the developers:

- fcomanuelgordilloc@gmail.com
- rosojeceb@gmail.com
- mrtprz@gmail.com

Hope you like it!!

## References

BOWTIE2:
- Langmead B. et al(2018). Scaling read aligners to hundreds of threads on general-purpose processors. Bioinformatics.
- Langmead B. et al (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods. 4;9(4):357-9.

FASTQC:
- Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online].

SAMTOOLS:
- Li H. et al (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.
MACS2:
- Zhang Y. et al (2008). Model-based Analysis of ChIP-Seq (MACS). Genome Biology, 9:R137.

HOMER: 
- Heinz S., Benner C., Spann N., Bertolino E. et al. (2010) Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol. Cell. May 28;38(4):576-589. 

R:   
- R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

RStudio:  
- RStudio Team (2020). RStudio: Integrated Development Environment for R. RStudio, PBC, Boston, MA. URL http://www.rstudio.com/.

Bioconductor packages for R analysis:
- ChIPseeker: Yu G. et al (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics 2015, 31(14):2382-2383
- TxDb.Athaliana.BioMart.plantsmart28: Carlson M, Maintainer BP (2015). TxDb.Athaliana.BioMart.plantsmart28: Annotation package for TxDb object(s).
- Org.at.tair.db: Carlson M. (2019). org.At.tair.db: Genome wide annotation for Arabidopsis. R package version 3.8.2.
- DO.db: Li J (2015). DO.db: A set of annotation maps describing the entire Disease Ontology. 
- ClusterProfiler: Yu G. et al (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5):284-287
- Pathview: Luo, W and Brouwer C (2013). Pathview: an R/Bioconductor package for pathway-based data integration and visualization. Bioinformatics, 29(14), 1830-1831.

KEGG:
- Kanehisa, M. and Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30.
- Kanehisa, M (2019). Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951.
- Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M., and Tanabe, M. (2021). KEGG: integrating viruses and cellular organisms. Nucleic Acids Res. 49, D545-D551.
