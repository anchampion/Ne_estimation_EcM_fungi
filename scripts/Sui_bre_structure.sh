# In the following, abbreviations correspond to the following populations
# :

# cal = California
# nocal = Inland (all except California)
# mocal = Mountain California
# cocal = Coastal California
# wycol = Wyoming + Colorado
# minwhi = Minnesota + Canada (Whitecourt)

module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.16

# For the "all" dataset, we want to have a representative sample across locations (3 from Mocal, 3 from Cocal, 3 from Wycol, 2 from Minwhi, 1 from Canada (Castle Rock)) 
# --> in each site, individuals are randomly chosen using the sample() function in R.

# Create a text file with the individuals' IDs (one per line)
gedit individuals.txt

# Create an output file containing only the individuals choosen
vcftools --vcf all20.recode.vcf --keep individuals.txt --recode --recode-INFO-all --out all20.12
rm *.log

# Then, we need to sample randomly 12 individuals from the different groups that contain more than 12 individuals : cal, nocal and mocal.

# Take 12 individuals randomly
vcftools --vcf cal20.recode.vcf --max-indv 12 --recode --recode-INFO-all --out cal20.12
vcftools --vcf mocal20.recode.vcf --max-indv 12 --recode --recode-INFO-all --out mocal20.12
vcftools --vcf nocal20.recode.vcf --max-indv 12 --recode --recode-INFO-all --out nocal20.12
rm *.log

# The 2 other pops have already 12 or 11 individuals, so we'll just rename them
cp minwhi.recode.vcf minwhi20.12.recode.vcf
cp wycol.recode.vcf wycol20.12.recode.vcf

# We also want to explore what happens with smaller subsets of SNPs, let's try 2000 SNPs.

# Generate files in the "plink" format, with extension .ped and .map, but we'll use only the .map files
mkdir data
vcftools --vcf all20.12.recode.vcf --plink --out ./data/all
vcftools --vcf mocal20.12.recode.vcf  --plink --out ./data/mocal 
vcftools --vcf cal20.12.recode.vcf --plink --out ./data/cal
vcftools --vcf nocal20.12.recode.vcf  --plink --out ./data/nocal
vcftools --vcf minwhi20.12.recode.vcf  --plink --out ./data/minwhi
vcftools --vcf wycol20.12.recode.vcf  --plink --out ./data/wycol

# deleting the log files that we don't need
rm ./data/*log

# We use the map file to extract the list of "scaffolds" and their variant sites (SNPs)
# saving in the newly created directory "snplist"
mkdir snplist
cut -f 2 ./data/all.map > ./snplist/all.snps
cut -f 2 ./data/cal.map > ./snplist/cal.snps
cut -f 2 ./data/mocal.map > ./snplist/mocal.snps
cut -f 2 ./data/nocal.map > ./snplist/nocal.snps
cut -f 2 ./data/minwhi.map > ./snplist/minwhi.snps
cut -f 2 ./data/wycol.map > ./snplist/wycol.snps

# From these lists, we subsample 2000 SNPs randomly using "shuffle = shuf"
# For each of the individual groups, we create 10 smaller subsets of SNPs:

# all
for i in {1..10}; do
 shuf ./snplist/all.snps | head -n 2000 | sort > ./snplist/all.subset2K$i.snps       
done

# cal
for i in {1..10}; do
 shuf ./snplist/cal.snps | head -n 2000 | sort > ./snplist/cal.subset2K$i.snps       
done

# mocal
for i in {1..10}; do
 shuf ./snplist/mocal.snps | head -n 2000 | sort > ./snplist/mocal.subset2K$i.snps       
done

# nocal
for i in {1..10}; do
 shuf ./snplist/nocal.snps | head -n 2000 | sort > ./snplist/nocal.subset2K$i.snps       
done

# minwhi
for i in {1..10}; do
 shuf ./snplist/minwhi.snps | head -n 2000 | sort > ./snplist/minwhi.subset2K$i.snps       
done

# minwhi
for i in {1..10}; do
 shuf ./snplist/wycol.snps | head -n 2000 | sort > ./snplist/wycol.subset2K$i.snps       
done

# we need to replace the ":" with a tab

# all
for i in {1..10}; do
 tr ':' '\t' < ./snplist/all.subset2K$i.snps > ./snplist/all.subset2.$i.snps    
done

# cal
for i in {1..10}; do
 tr ':' '\t' < ./snplist/cal.subset2K$i.snps > ./snplist/cal.subset2.$i.snps    
done

# mocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist/mocal.subset2K$i.snps > ./snplist/mocal.subset2.$i.snps    
done

# nocal
for i in {1..10}; do
 tr ':' '\t' < ./snplist/nocal.subset2K$i.snps > ./snplist/nocal.subset2.$i.snps    
done

# minwhi
for i in {1..10}; do
 tr ':' '\t' < ./snplist/minwhi.subset2K$i.snps > ./snplist/minwhi.subset2.$i.snps    
done

# minwhi
for i in {1..10}; do
 tr ':' '\t' < ./snplist/wycol.subset2K$i.snps > ./snplist/wycol.subset2.$i.snps    
done

# notice the last files do not have "K" in their name. 
# and then rename them

for i in {1..10}; do
 mv ./snplist/all.subset2.$i.snps ./snplist/all.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist/cal.subset2.$i.snps ./snplist/cal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist/mocal.subset2.$i.snps ./snplist/mocal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist/nocal.subset2.$i.snps ./snplist/nocal.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist/minwhi.subset2.$i.snps ./snplist/minwhi.subset2K$i.snps   
done

for i in {1..10}; do
 mv ./snplist/wycol.subset2.$i.snps ./snplist/wycol.subset2K$i.snps   
done

# with vcftools, we create the vcf files with the new subset of SNPs

for i in {1..10}; do
vcftools --vcf all20.12.recode.vcf --positions ./snplist/all.subset2K$i.snps --recode --out ./data/all2K.$i
done

for i in {1..10}; do
vcftools --vcf cal20.12.recode.vcf --positions ./snplist/cal.subset2K$i.snps --recode --out ./data/cal2K.$i
done

for i in {1..10}; do
vcftools --vcf mocal20.12.recode.vcf --positions ./snplist/mocal.subset2K$i.snps --recode --out ./data/mocal2K.$i
done

for i in {1..10}; do
vcftools --vcf nocal20.12.recode.vcf --positions ./snplist/nocal.subset2K$i.snps --recode --out ./data/nocal2K.$i
done

for i in {1..10}; do
vcftools --vcf minwhi20.12.recode.vcf --positions ./snplist/minwhi.subset2K$i.snps --recode --out ./data/minwhi2K.$i
done

for i in {1..10}; do
vcftools --vcf wycol20.12.recode.vcf --positions ./snplist/wycol.subset2K$i.snps --recode --out ./data/wycol2K.$i
done