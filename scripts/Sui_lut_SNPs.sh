module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.16

# checking missing data first:
vcftools --vcf gvcf.filtered.2020-05-23dik.recode.vcf --remove-indels --max-missing 0.8 --recode --out Sui_lut

# we generate files in the "plink" format, with extension .ped and .map, but we'll use only the .map files

vcftools --vcf Sui_lut.recode.vcf --plink --out Sui_lut 

# deleting the the log files that we don't need

rm *log

# We use the map file to extract the list of "scaffolds" and their variant sites (SNPs)
# saving in the newly created directory "snplist"
mkdir snplist
cut -f 2 Sui_lut.map > ./snplist/Sui_lut.snps

# From these lists, we subsample 2000 SNPs randomly using "shuffle = shuf"
# For each of the individual groups, we create 10 smaller subsets of SNPs:

for i in {1..10}; do
 shuf ./snplist/Sui_lut.snps | head -n 2000 | sort > ./snplist/Sui_lut.subset2K$i.snps       
done

# we need to replace the ":" with a tab

for i in {1..10}; do
 tr ':' '\t' < ./snplist/Sui_lut.subset2K$i.snps > ./snplist/Sui_lut.subset2.$i.snps    
done

# notice the last files do not have "K" in their name. 
# and then rename them

for i in {1..10}; do
 mv ./snplist/Sui_lut.subset2.$i.snps ./snplist/Sui_lut.subset2K$i.snps   
done

# with vcftools, we create the vcf files with the new subset of SNPs

for i in {1..10}; do
vcftools --vcf Sui_lut.recode.vcf --positions ./snplist/Sui_lut.subset2K$i.snps --recode --out ./data/Sui_lut2K.$i
done

## Let's do the same for 20K SNPs subsets

# From these lists, we subsample 20000 SNPs randomly using "shuffle = shuf"
# For each of the individual groups, we create 10 smaller subsets of SNPs:

for i in {1..10}; do
 shuf ./snplist/Sui_lut.snps | head -n 20000 | sort > ./snplist/Sui_lut.subset20K$i.snps       
done

# we need to replace the ":" with a tab

for i in {1..10}; do
 tr ':' '\t' < ./snplist/Sui_lut.subset20K$i.snps > ./snplist/Sui_lut.subset20.$i.snps    
done

# notice the last files do not have "K" in their name. 
# and then rename them 

for i in {1..10}; do
 mv ./snplist/Sui_lut.subset20.$i.snps ./snplist/Sui_lut.subset20K$i.snps   
done

# with vcftools, we create the vcf files with the new subset of SNPs

for i in {1..10}; do
vcftools --vcf Sui_lut.recode.vcf --positions ./snplist/Sui_lut.subset20K$i.snps --recode --out ./data/Sui_lut20K.$i
done