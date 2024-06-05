module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.16

# Sample 10 replicates of 2 000 SNPs in each pop in order to compare estimations together.
# We have to start with the original dataset ("XXX.recode.vcf") and sampling 12 individuals in each population.

# Create an output file containing only the individuals choosen in "all"
vcftools --vcf all.recode.vcf --keep individuals.txt --recode --recode-INFO-all --out all.12
rm *.log

# Then, we need to sample randomly 12 individuals from the different groups that contain more than 12 individuals : cal, nocal and mocal.

# Take 12 individuals randomly
vcftools --vcf cal.recode.vcf --max-indv 12 --recode --recode-INFO-all --out cal.12
vcftools --vcf mocal.recode.vcf --max-indv 12 --recode --recode-INFO-all --out mocal.12
vcftools --vcf nocal.recode.vcf --max-indv 12 --recode --recode-INFO-all --out nocal.12
rm *.log

# The 2 other pops have already 12 or 11 individuals, so we'll just rename them
cp minwhi.recode.vcf minwhi.12.recode.vcf
cp wycol.recode.vcf wycol.12.recode.vcf

# we first generate files in the "plink" format, with extension .ped and .map, but we'll use only the .map files

vcftools --vcf all.12.recode.vcf --plink --out ./data/all 
vcftools --vcf mocal.12.recode.vcf  --plink --out ./data/mocal 
vcftools --vcf cal.12.recode.vcf --plink --out ./data/cal 
vcftools --vcf minwhi.12.recode.vcf --plink --out ./data/minwhi
vcftools --vcf wycol.12.recode.vcf  --plink --out ./data/wycol 
vcftools --vcf nocal.12.recode.vcf  --plink --out ./data/nocal 

# deleting the log files that we don't need
rm ./data/*log

# We use the map file to extract the list of "scaffolds" and their variant sites (SNPs)
# saving in the newly created directory "snplist2"
mkdir snplist2
cut -f 2 ./data/all.map > ./snplist2/all.snps
cut -f 2 ./data/cal.map > ./snplist2/cal.snps
cut -f 2 ./data/mocal.map > ./snplist2/mocal.snps
cut -f 2 ./data/nocal.map > ./snplist2/nocal.snps
cut -f 2 ./data/minwhi.map > ./snplist2/minwhi.snps
cut -f 2 ./data/wycol.map > ./snplist2/wycol.snps

# From these lists, we subsample 2000 SNPs randomly using "shuffle = shuf"
# For each of the individual groups, we create 10 smaller subsets of SNPs:

# all
for i in {1..10}; do
 shuf ./snplist2/all.snps | head -n 2000 | sort > ./snplist2/all.subset2K$i.snps       
done

# cal
for i in {1..10}; do
 shuf ./snplist2/cal.snps | head -n 2000 | sort > ./snplist2/cal.subset2K$i.snps       
done

# mocal
for i in {1..10}; do
 shuf ./snplist2/mocal.snps | head -n 2000 | sort > ./snplist2/mocal.subset2K$i.snps       
done

# nocal
for i in {1..10}; do
 shuf ./snplist2/nocal.snps | head -n 2000 | sort > ./snplist2/nocal.subset2K$i.snps       
done

# minwhi
for i in {1..10}; do
 shuf ./snplist2/minwhi.snps | head -n 2000 | sort > ./snplist2/minwhi.subset2K$i.snps       
done

# mycol
for i in {1..10}; do
 shuf ./snplist2/wycol.snps | head -n 2000 | sort > ./snplist2/wycol.subset2K$i.snps       
done

# we need to replace the ":" with a tab

# all
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/all.subset2K$i.snps > ./snplist2/all.subset2.$i.snps    
done

# cal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/cal.subset2K$i.snps > ./snplist2/cal.subset2.$i.snps    
done

# mocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/mocal.subset2K$i.snps > ./snplist2/mocal.subset2.$i.snps    
done

# nocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/nocal.subset2K$i.snps > ./snplist2/nocal.subset2.$i.snps    
done

# minwhi
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/minwhi.subset2K$i.snps > ./snplist2/minwhi.subset2.$i.snps    
done

# wycol
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/wycol.subset2K$i.snps > ./snplist2/wycol.subset2.$i.snps    
done

# notice the last files do not have "K" in their name. 
# and then rename them

for i in {1..10}; do
 mv ./snplist2/all.subset2.$i.snps ./snplist2/all.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/cal.subset2.$i.snps ./snplist2/cal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/mocal.subset2.$i.snps ./snplist2/mocal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/nocal.subset2.$i.snps ./snplist2/nocal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/minwhi.subset2.$i.snps ./snplist2/minwhi.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/wycol.subset2.$i.snps ./snplist2/wycol.subset2K$i.snps   
done

# with vcftools, we create the vcf files with the new subset of SNPs

for i in {1..10}; do
vcftools --vcf all.12.recode.vcf --positions ./snplist2/all.subset2K$i.snps --recode --out ./data/all2K.$i
done

for i in {1..10}; do
vcftools --vcf cal.12.recode.vcf --positions ./snplist2/cal.subset2K$i.snps --recode --out ./data/cal2K.$i
done

for i in {1..10}; do
vcftools --vcf mocal.12.recode.vcf --positions ./snplist2/mocal.subset2K$i.snps --recode --out ./data/mocal2K.$i
done

for i in {1..10}; do
vcftools --vcf nocal.12.recode.vcf --positions ./snplist2/nocal.subset2K$i.snps --recode --out ./data/nocal2K.$i
done

for i in {1..10}; do
vcftools --vcf minwhi.12.recode.vcf --positions ./snplist2/minwhi.subset2K$i.snps --recode --out ./data/minwhi2K.$i
done

for i in {1..10}; do
vcftools --vcf wycol.12.recode.vcf --positions ./snplist2/wycol.subset2K$i.snps --recode --out ./data/wycol2K.$i
done

# From the same lists, we subsample 20000 SNPs randomly using "shuffle = shuf"
# For each of the individual groups, we create 10 smaller subsets of SNPs:

# all
for i in {1..10}; do
 shuf ./snplist2/all.snps | head -n 20000 | sort > ./snplist2/all.subset20K$i.snps       
done

# cal
for i in {1..10}; do
 shuf ./snplist2/cal.snps | head -n 20000 | sort > ./snplist2/cal.subset20K$i.snps       
done

# mocal
for i in {1..10}; do
 shuf ./snplist2/mocal.snps | head -n 20000 | sort > ./snplist2/mocal.subset20K$i.snps       
done

# nocal
for i in {1..10}; do
 shuf ./snplist2/nocal.snps | head -n 20000 | sort > ./snplist2/nocal.subset20K$i.snps       
done

# minwhi
for i in {1..10}; do
 shuf ./snplist2/minwhi.snps | head -n 20000 | sort > ./snplist2/minwhi.subset20K$i.snps       
done

# mycol
for i in {1..10}; do
 shuf ./snplist2/wycol.snps | head -n 20000 | sort > ./snplist2/wycol.subset20K$i.snps       
done

# we need to replace the ":" with a tab

# all
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/all.subset20K$i.snps > ./snplist2/all.subset20.$i.snps    
done

# cal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/cal.subset20K$i.snps > ./snplist2/cal.subset20.$i.snps    
done

# mocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/mocal.subset20K$i.snps > ./snplist2/mocal.subset20.$i.snps    
done

# nocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/nocal.subset20K$i.snps > ./snplist2/nocal.subset20.$i.snps    
done

# minwhi
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/minwhi.subset20K$i.snps > ./snplist2/minwhi.subset20.$i.snps    
done

# wycol
for i in {1..10}; do
 tr ':' '\t' < ./snplist2/wycol.subset20K$i.snps > ./snplist2/wycol.subset20.$i.snps    
done

# notice the last files do not have "K" in their name. 
# and then rename them 

for i in {1..10}; do
 mv ./snplist2/all.subset20.$i.snps ./snplist2/all.subset20K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/cal.subset20.$i.snps ./snplist2/cal.subset20K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/mocal.subset20.$i.snps ./snplist2/mocal.subset20K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/nocal.subset20.$i.snps ./snplist2/nocal.subset20K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/minwhi.subset20.$i.snps ./snplist2/minwhi.subset20K$i.snps   
done

for i in {1..10}; do
 mv ./snplist2/wycol.subset20.$i.snps ./snplist2/wycol.subset20K$i.snps   
done

# with vcftools, we create the vcf files with the new subset of SNPs

for i in {1..10}; do
vcftools --vcf all.12.recode.vcf --positions ./snplist2/all.subset20K$i.snps --recode --out ./data/all20K.$i
done

for i in {1..10}; do
vcftools --vcf cal.12.recode.vcf --positions ./snplist2/cal.subset20K$i.snps --recode --out ./data/cal20K.$i
done

for i in {1..10}; do
vcftools --vcf mocal.12.recode.vcf --positions ./snplist2/mocal.subset20K$i.snps --recode --out ./data/mocal20K.$i
done

for i in {1..10}; do
vcftools --vcf nocal.12.recode.vcf --positions ./snplist2/nocal.subset20K$i.snps --recode --out ./data/nocal20K.$i
done

for i in {1..10}; do
vcftools --vcf minwhi.12.recode.vcf --positions ./snplist2/minwhi.subset20K$i.snps --recode --out ./data/minwhi20K.$i
done

for i in {1..10}; do
vcftools --vcf wycol.12.recode.vcf --positions ./snplist2/wycol.subset20K$i.snps --recode --out ./data/wycol20K.$i
done