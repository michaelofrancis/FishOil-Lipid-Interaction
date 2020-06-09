###Michael Francis 06-09-2020
###This script runs a GWAS interaction model on four blood lipid phenotypes using fish oil as an interaction term.
###http://toby.freeshell.org/software/quicktest.shtml


###Set Parameters ==============================

#Today's date
now=$(date +"%m_%d_%Y")

#Set phenotypes
phenotypes=("LDL" "HDL" "Tot_Chol" "TAGs")

#Set imputed genotype data input directory
genoindir=("")

#Set directory in which to generate scripts (one subfolder will be created).
#Ideally this is the directory in which this master script is located.
scriptdir=("")

#Set output directory (subfolders will be created)
outdir=("")

#Set Phenotype Dir (Absolute path)
phenodir=("")

#Set queue's to run on
mainq=("")

chr=(22 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)

#End Set Parameters====================================

#~~~~~~~~~~~~~~~~~~~~~~~BEGIN LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for j in ${phenotypes[@]}
        do
	for i in ${chr[*]}
		do

mkdir -p "$scriptdir"/"$j"
cd "$scriptdir"/"$j"
mkdir -p "$outdir"/"$j"

##NOTE: i's are initially made as question marks then replaced with a sed command after.

#Generate scripts=-=-=-=-=-=-=-=-=--=-=-=-
echo "#PBS -S /bin/bash
#PBS -N quicktest-hm200-"$j"-chr"$i"
#PBS -q "$mainq"
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00:00
#PBS -l mem=200gb

module load R/3.6.2-foss-2018a-X11-20180131-GACRC
export PATH=\${HOME}/bin/quicktest/1.2:\${PATH}

cd "$scriptdir"/"$j"
mkdir -p "$outdir"/"$j"/

quicktest \
--pheno "$phenodir"/"$j"-phenochr"$i"sorted-test-03252020.txt2.txt \
--snptest \
--geno "$genoindir"/"$j"/"$j".chr"$i".gen.gz \
--npheno "$j" \
--ncovar FishOil2 \
--ncovar Sex --ncovar Age --ncovar Townsend --ncovar weekly_oily_fish --ncovar BMI \
--ncovar PCA1 --ncovar PCA2 --ncovar PCA3 --ncovar PCA4 --ncovar PCA5 \
--ncovar PCA6 --ncovar PCA7 --ncovar PCA8 --ncovar PCA9 --ncovar PCA10 \
--method-interaction \
--method-robust \
--out "$outdir"/"$j"/chr"$i".hm200-04122020.out

" > "$scriptdir"/"$j"/"$j"_chr"$i"_Quicktest_UKB.sh

qsub "$scriptdir"/"$j"/"$j"_chr"$i"_Quicktest_UKB.sh

done
done
