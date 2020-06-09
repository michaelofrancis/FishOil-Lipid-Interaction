###Michael Francis, 06-09-2020
###This script runs a meta-analysis on QuickTest interaction GWAS output files. 
###METAL is patched with interaction functions from A. Manning
###https://genome.sph.umich.edu/wiki/Meta_Analysis_of_SNPxEnvironment_Interaction

#STAGE 1+2 META-ANALYSIS

#VERBOSE ON

SCHEME INTERACTION

#Describe and process files
SEPARATOR WHITESPACE
MARKER   SNP
ALLELE   alleleA alleleB
EFFECT   snp.beta
STDERR   snp.se
INTEFFECT interaction.beta
INTSTDERR interaction.se
INTCOV	cov.snp.interaction
WEIGHT N
FREQ ALT_FREQS
MINMAXFREQ ON
AVERAGEFREQ ON
GENOMICCONTROL ON

PROCESS #path to UK Biobank discovery significant variants QuickTest output file(s) (1*10^-6)  
PROCESS #ARIC QuickTest output file(s)


OUTFILE LDL.stage3 .txt

# Execute meta-analysis
ANALYZE
