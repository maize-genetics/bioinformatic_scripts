# Notebooks for the Practical Haplotype Graph

The notebooks included here pertain to building a new practical haplotype graph. Use build-a-phg.ipynb to walk through the steps required to create a new PHG from scratch. As you generate GVCF and SAM files, use the respective Metrics notebooks to extract information such as read coverage, number of SNPs, and indel size distribution. 

The build-a-phg walkthrough assumes that you have installed docker, faSize, and kotlinc, and that all are on your path. It also assumes that you have a reference .fasta file and associated gene annotations in .gff format, and that for each assembly you want to use to generate haplotypes you have a .maf file from aligning the assembly to your reference.

The metric notebooks require .gvcf and .sam files, respectively, and they also require that you have installed docker. Additionally, they use a number of python packages related to data manipulation and graphing. Many of these are common packages that you may already have, but check the header of each notebook before you begin for the complete list of packages to install.  