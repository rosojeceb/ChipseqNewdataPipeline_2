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
	- clusterProfiler

To run the program, please make sure you filla ll the test_params.txt parameters correctly and then launch both with the ./ command. Whenever you launch ChipseqNewdataPipeline without parameters, it will also ask for this information.

If you have any further questions, do not hesitate to ask the developers:

	- fcomanuelgordilloc@gmail.com
	- rosojeceb@gmail.com
	- mrtprz@gmail.com

Hope you like it!!

 References

