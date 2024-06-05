# Ne_estimation_EcM_fungi
 Code and Workflow - M2 Internship

Analyses for Ne estimation in ectomycorrhizal fungi
================
Anouck Champion
2024-06-05

This document describes the analyses carried out in R and Genotoul
bioinfo to determine the effective population size (Ne) in 4 species of
ectomycorrhizal fungi : Boletus edulis, Suillus brevipes, Suillus luteus
and Tuber melanosporum. These analyses include a filtering phase of SNPs
and microsatellite data, and an analysis part dedicated to the
identification of bias in Ne estimation.

## Packages loading

## Effect of ploidy on Ne estimations

Dataset used : *Tuber melanosporum* from Taschen et al. 2016

We want to compare Ne estimations (and other metrics : He, Fis) for 3
different approaches :

- 1 : Use the “zygotes” (diploid data) provided in the article
- 2 : Duplicate haploid data to create artificial diploid genotypes
  (homozygotes only)
- 3 : Combine haploid genotypes to create artificial diploid genotypes

### Filtering step

Before accounting for the impact of ploidy, we must check 2 things :

- are they missing data in the dataset ?
- are they clones in the dataset ?
- is the population homogenous or structured ?

#### Diploid data (zygotes)

First exploration of the data, to see if we can combine several
locations in one population.

``` r
data<-read.genepop(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "zygotes", 
                        "zygotes_all_noclones.gen"), ncode = 3L, quiet = TRUE)

# Clones were already removed from this dataset, so we move to structure
# Assessing structure using DAPC from the adegenet package
data
data@pop
dapc0 <- dapc(data, n.pca = 5)
scatter(dapc0)
# We will use only the SB1 and SB2 sites, and consider them as 1 homogenous population, with regard to DAPC
```

Then, starting from the original dataset. This is an example for SB1,
but the same code was used for SB2.

``` r
# SB1
zyg_SB1<-read.genepop(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "zygotes", 
                        "zygotes_SB1.gen"), ncode = 3L, quiet = TRUE)
zyg_SB1
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 20 individuals; 13 loci; 35 alleles; size: 17.6 Kb
    ## 
    ##  // Basic content
    ##    @tab:  20 x 35 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-5)
    ##    @loc.fac: locus factor for the 35 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: read.genepop(file = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", 
    ##     "zygotes", "zygotes_SB1.gen"), ncode = 3L, quiet = TRUE)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 20-20)

``` r
# Missing data ? (with poppr package)
missing.loc<-missingno(zyg_SB1, type = "loci", cutoff = 0.20);missing.loc
```

    ## 
    ## Found 53 missing values.
    ## 
    ## 2 loci contained missing values greater than 20%
    ## 
    ## Removing 2 loci: Tm2, Tm98

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 20 individuals; 11 loci; 29 alleles; size: 15.2 Kb
    ## 
    ##  // Basic content
    ##    @tab:  20 x 29 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-5)
    ##    @loc.fac: locus factor for the 29 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: .local(x = x, i = i, j = j, drop = drop)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 20-20)

``` r
missing.ind<-missingno(missing.loc, type = "genotype", cutoff = 0.20);missing.ind
```

    ## 
    ##  No genotypes with missing values above 20% found.

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 20 individuals; 11 loci; 29 alleles; size: 15.2 Kb
    ## 
    ##  // Basic content
    ##    @tab:  20 x 29 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-5)
    ##    @loc.fac: locus factor for the 29 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: .local(x = x, i = i, j = j, drop = drop)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 20-20)

``` r
zyg_SB1<-missing.ind;zyg_SB1
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 20 individuals; 11 loci; 29 alleles; size: 15.2 Kb
    ## 
    ##  // Basic content
    ##    @tab:  20 x 29 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-5)
    ##    @loc.fac: locus factor for the 29 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: .local(x = x, i = i, j = j, drop = drop)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 20-20)

``` r
# Clones ? (with Rclone method)

# Reformat the data for Rclone
data.df <- genind2df(zyg_SB1)
popvec <- data.df[,1] # pop info (if present)
dataset <- data.df[,2:ncol(data.df)]# dataset
dataset <- as.data.frame(apply(dataset, 2, function(x) replace(x, is.na(x), "000000"))) # replace NAs by 000000
data2 <- convert_GC(dataset, 3) # we use "3" because this is the length of the allele
row.names(data2) <- seq_len(nrow(data2)) # make rownames consecutive integers starting from 1 to avoid bugs
head(data2)
```

    ##   me02_1 me02_2 me11_1 me11_2 me13_1 me13_2 me14_1 me14_2 Tm1_1 Tm1_2 Tm9_1
    ## 1    155    155    302    302    118    118    142    142   342   342   326
    ## 2    155    155    302    302    118    118    142    142   342   342   326
    ## 3    155    155    302    302    118    118    136    142   342   342   326
    ## 4    155    155    302    302    118    118    142    142   342   342   326
    ## 5    155    155    302    302    118    118    142    142   342   342   326
    ## 6    155    155    302    302    118    118    142    142   000   000   326
    ##   Tm9_2 Tm22_1 Tm22_2 Tm21_1 Tm21_2 Tm75_1 Tm75_2 Tm269_1 Tm269_2 Tm127_1
    ## 1   326    322    322    313    313    342    342     000     000     182
    ## 2   326    322    328    313    316    342    342     371     371     182
    ## 3   326    322    322    316    316    342    342     371     371     182
    ## 4   326    322    328    316    316    342    342     371     371     182
    ## 5   326    322    322    313    313    342    342     000     000     182
    ## 6   326    322    322    000    000    342    342     371     371     182
    ##   Tm127_2
    ## 1     182
    ## 2     182
    ## 3     192
    ## 4     192
    ## 5     182
    ## 6     182

``` r
# We can write this dataframe to a new file
write.table(data2, 
            file = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "zygotes","SB1_Rclone.txt"),
            row.names = FALSE,
            sep = "\t")

# Discrimination of MLG
# List unique alleles per locus
list_all_tab(data2)
```

    ##   locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9
    ## 1     155     302     118     142     342     326     322     313     342
    ## 2     161     296             136     000     301     328     316     326
    ## 3     171     284             148     343     330     335     000        
    ## 4     165                                                     289        
    ## 5                                                             256        
    ## 6                                                             310        
    ##   locus_10 locus_11
    ## 1      000      182
    ## 2      371      192
    ## 3                  
    ## 4                  
    ## 5                  
    ## 6

``` r
# List MLG
MLG_tab(data2)
```

    ##    unit_1 unit_2 unit_3
    ## 1       1      5       
    ## 2       2              
    ## 3       3              
    ## 4       4              
    ## 5       6              
    ## 6       7     10     13
    ## 7       8              
    ## 8       9     17     19
    ## 9      11              
    ## 10     12              
    ## 11     14              
    ## 12     15              
    ## 13     16              
    ## 14     18              
    ## 15     20

``` r
MLGlist <- MLG_list(data2)
# Allelic frequencies :
freq_RR(data2)
```

    ##       locus allele  freq
    ## 1   locus_1    155 0.925
    ## 2   locus_1    161 0.025
    ## 3   locus_1    165 0.025
    ## 4   locus_1    171 0.025
    ## 5   locus_2    284 0.025
    ## 6   locus_2    296 0.050
    ## 7   locus_2    302 0.925
    ## 8   locus_3    118 1.000
    ## 9   locus_4    136 0.025
    ## 10  locus_4    142 0.950
    ## 11  locus_4    148 0.025
    ## 12  locus_5    000 0.050
    ## 13  locus_5    342 0.900
    ## 14  locus_5    343 0.050
    ## 15  locus_6    301 0.025
    ## 16  locus_6    326 0.950
    ## 17  locus_6    330 0.025
    ## 18  locus_7    322 0.700
    ## 19  locus_7    328 0.275
    ## 20  locus_7    335 0.025
    ## 21  locus_8    000 0.050
    ## 22  locus_8    256 0.025
    ## 23  locus_8    289 0.025
    ## 24  locus_8    310 0.025
    ## 25  locus_8    313 0.725
    ## 26  locus_8    316 0.150
    ## 27  locus_9    326 0.050
    ## 28  locus_9    342 0.950
    ## 29 locus_10    000 0.200
    ## 30 locus_10    371 0.800
    ## 31 locus_11    182 0.950
    ## 32 locus_11    192 0.050

``` r
# OR using poppr
# Calculating genotypic diversity
poppr(zyg_SB1) # 20 ramets but only 15 genets (MLG)
```

    ##    Pop  N MLG eMLG SE   H    G lambda   E.5  Hexp    Ia  rbarD    File
    ## 1 260P 20  15   15  0 2.6 11.8  0.915 0.867 0.148 0.324 0.0457 zyg_SB1

``` r
# Plot the MLG per sampling site
P.tab <- mlg.table(zyg_SB1)
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Determination of MLL

#genetic distances computation, distance on allele differences:
respop <- genet_dist(data2);respop # max 10 differences between genotypes
```

    ## $distance_matrix
    ##     1  2  3  4  5  6  7  8  9 10 11 12 13 14
    ## 2   4                                       
    ## 3   6  4                                    
    ## 4   6  2  2                                 
    ## 5   6  5  6  6                              
    ## 6   4  2  6  4  6                           
    ## 7   4  2  3  3  4  4                        
    ## 8   2  2  4  4  4  2  2                     
    ## 9   5  6 10  8 10  5  8  7                  
    ## 10  3  7  8  9  9  7  7  5  6               
    ## 11  7  5  9  7  9  4  7  5  8  9            
    ## 12  3  2  4  4  4  3  2  1  7  6  6         
    ## 13  3  3  5  5  4  3  3  1  8  6  6  2      
    ## 14  5  5  7  7  6  5  5  3 10  8  8  4  2   
    ## 15  3  3  5  5  5  3  3  1  8  6  6  2  2  3

``` r
ressim <- genet_dist_sim(data2, nbrepeat = 1000) #theoretical distribution : sexual reproduction
```

    ## [1] "Number of MLG sim = 233"

``` r
ressimWS <- genet_dist_sim(data2, genet = TRUE, nbrepeat = 1000) #idem, without selfing
```

    ## [1] "Number of MLG sim = 291"

``` r
#graph prep.:
p1 <- hist(respop$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), main = "popsim", 
           xlab = "Genetic distances", breaks = seq(0, max(respop$distance_matrix)+1, 1))
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
p2 <- hist(ressim$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), main = "popSR",
           xlab = "Genetic distances", breaks = seq(0, max(ressim$distance_matrix)+1, 1))
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
p3 <- hist(ressimWS$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3),
           main = "popSRWS", xlab = "Genetic distances", breaks = seq(0, max(ressimWS$distance_matrix)+1, 1))
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
limx <- max(max(respop$distance_matrix), max(ressim$distance_matrix), max(ressimWS$distance_matrix))

#graph superposition:
plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), main = "", xlab = "Genetic distances")
plot(p2, col = rgb(0.7,0.9,1,0.5), freq = FALSE, add = TRUE)
plot(p3, col = rgb(0.9,0.5,1,0.3), freq = FALSE, add = TRUE)

#adding a legend:
leg.txt <- c("original data","simulated data", "without selfing")
col <- c(rgb(0,0.4,1,1), rgb(0.7,0.9,1,0.5), rgb(0.9,0.5,1,0.3))
legend("top", fill = col, leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white")
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

``` r
#determining alpha2 (= highest genetic distance (i.e. max number of differences) between genotypes that we condider for MLL)
# Can be determined graphically (cutoff) or by looking directly at the data :
table(respop$distance_matrix)
```

    ## 
    ##  1  2  3  4  5  6  7  8  9 10 
    ##  3 14 15 16 15 16 10  8  5  3

``` r
# No MLL here --> keep MLG (i.e. alpha2 = 0)

MLLlist <- MLL_generator(data2, alpha2 = 0) # Change alpha2 value
#or (if problem)
# res <- genet_dist(data2, alpha2 = 2)
# MLLlist <- MLL_generator2(potential_clones = res$potential_clones, res_mlg = MLGlist)

MLLlist
```

    ## [[1]]
    ## [1] 1 5
    ## 
    ## [[2]]
    ## [1] 2
    ## 
    ## [[3]]
    ## [1] 3
    ## 
    ## [[4]]
    ## [1] 4
    ## 
    ## [[5]]
    ## [1] 6
    ## 
    ## [[6]]
    ## [1]  7 10 13
    ## 
    ## [[7]]
    ## [1] 8
    ## 
    ## [[8]]
    ## [1]  9 17 19
    ## 
    ## [[9]]
    ## [1] 11
    ## 
    ## [[10]]
    ## [1] 12
    ## 
    ## [[11]]
    ## [1] 14
    ## 
    ## [[12]]
    ## [1] 15
    ## 
    ## [[13]]
    ## [1] 16
    ## 
    ## [[14]]
    ## [1] 18
    ## 
    ## [[15]]
    ## [1] 20

``` r
# Extract one individual from each MLL 
new_ind <- sapply(MLLlist, `[[`, 1)
new_ind <- as.numeric(new_ind);new_ind
```

    ##  [1]  1  2  3  4  6  7  8  9 11 12 14 15 16 18 20

``` r
rownames(zyg_SB1@tab)<-c(1:20) # if ids are not set to 1,2,3,...n
selected_indices <- match(new_ind, rownames(zyg_SB1@tab));selected_indices
```

    ##  [1]  1  2  3  4  6  7  8  9 11 12 14 15 16 18 20

``` r
# Create a new genind object without the clonal individuals
data_noclones <- zyg_SB1[selected_indices, ]
# Check if the genind object is correct
data_noclones
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 15 individuals; 11 loci; 29 alleles; size: 14.2 Kb
    ## 
    ##  // Basic content
    ##    @tab:  15 x 29 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 1-5)
    ##    @loc.fac: locus factor for the 29 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: .local(x = x, i = i, j = j, drop = drop)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 15-15)

``` r
# When it's ok --> make this new genind the main genind object
zyg_SB1 <- data_noclones

# Structure ?

# Population structure has been already explored before. SB1 is a homogenous population.

# Export the final dataset that will be used after
genind_to_genepop(zyg_SB1, output = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", "zyg_SB1_correct.txt"))
```

    ## There is only one population in your dataset

``` r
# ==> Use the same code for SB2
# As only 11 loci were kept for SB1 and for the haploid dataset, we also keep the same 11 loci for SB2
locNames(zyg_SB2)
# We want to remove loci Tm2 and Tm98
toRemove <- c(7,12)
zyg_SB2 <- zyg_SB2[loc=-toRemove]
zyg_SB2
locNames(zyg_SB2)

genind_to_genepop(zyg_SB2, output = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", "zyg_SB2_correct.txt"))
```

#### Duplicated data (from haploids)

The same workflow is used here on the haploid dataset. Once the dataset
is filtered, the duplication is processed using the GenAlEx programme
implemented in Excel.

``` r
# SB1
haplo_SB1 <- read.genalex(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "haploid_data", 
                          "SB1.csv"), ploidy = 1, genclone = FALSE, sep = ";")
haplo_SB1

# Missing data ?
missing.loc<-missingno(haplo_SB1, type = "loci", cutoff = 0.20);missing.loc
missing.ind<-missingno(missing.loc, type = "genotype", cutoff = 0.20);missing.ind

haplo_SB1<-missing.ind;haplo_SB1

# Clones ? (with Rclone method - for haploid data)

data.df <- genind2df(haplo_SB1);View(data.df)
vechaplo <- data.df[,1] # pop info (if present)
dataset <- data.df[,2:ncol(data.df)]# dataset
haplodata <- dataset 
head(haplodata)

MLG_tab(haplodata)

respop <- genet_dist(haplodata, haploid = TRUE);respop
ressim <- genet_dist_sim(haplodata, haploid = TRUE, nbrepeat = 1000)

p1 <- hist(respop$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), main = "popsim", xlab = "Genetic distances", breaks = seq(0, max(respop$distance_matrix)+1, 1))
p2 <- hist(ressim$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), main = "popSR", xlab = "Genetic distances", breaks = seq(0, max(ressim$distance_matrix)+1, 1))
limx <- max(max(respop$distance_matrix), max(ressim$distance_matrix))

plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), main = "", xlab = "Genetic distances")
plot(p2, col = rgb(0.7,0.9,1,0.5), freq = FALSE, add = TRUE)

MLLlist <- MLL_generator(haplodata, haploid = TRUE, alpha2 = 0)
MLLlist

# Extract one individual from each MLL 
new_ind <- sapply(MLLlist, `[[`, 1)
new_ind <- as.numeric(new_ind);new_ind

rownames(haplo_SB1@tab)<-c(1:135) # if ids are not set to 1,2,3,...n
selected_indices <- match(new_ind, rownames(haplo_SB1@tab));selected_indices

# Create a new genind object without the clonal individuals
data_noclones <- haplo_SB1[selected_indices, ]
# Check if the genind object is correct
haplo_SB1 <- data_noclones

# Export to csv
genind2genalex(data_noclones, filename = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", 
                                              "haploid_data", "SB2_haplo_correct.csv"), sep = ";")

# After duplication via genalex, here is the final dataset
dup_SB1 <- read.genalex(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                             "dup_SB1_correct.csv"), ploidy = 2, genclone = FALSE, sep = ";")
dup_SB1

# After duplication via genalex, here is the final dataset
dup_SB1 <- read.genalex(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                             "dup_SB1_correct.csv"), ploidy = 2, genclone = FALSE, sep = ";")
dup_SB1

# ==> Use the same code for SB2
```

#### Combined data (from haploids)

As described in the methods, the second method used to deal with haploid
data consisted of pairing randomly haploid genotypes to create diploids.
Warning : this method is called “combined” or “combination” in this
scripts but refers to the “pairing” method described in the methods.
Here the same haploid datasets (SB1 & SB2) as for the “duplication
method” are used.

``` r
# SB1
haplo_SB1
# Sample at random 18 (or 12) individuals from the 19 (we need an even number !)
# Number of individuals in the object
n_ind <- nrow(haplo_SB1@tab)

# Set the seed for reproducibility
set.seed(123)

# Sample 18 random individuals without replacement
sampled_indices <- sample(n_ind, size = 18, replace = FALSE)

# Subset the genind object based on the sampled indices
sampled_genind <- haplo_SB1[sampled_indices, ]

# Export to csv
genind2genalex(sampled_genind, filename = here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", 
                                              "haploid_data", "subsamples_for_combination", "SB1_haplo_correct_n18.csv"), sep = ";")

# After combination via genalex, here is the final dataset
comb_SB1 <- read.genalex(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                             "comb_SB1_correct.csv"), ploidy = 2, genclone = FALSE, sep = ";")

comb_SB1

# After combination via genalex, here is the final dataset
comb_SB1 <- read.genalex(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                              "comb_SB2_correct.csv"), ploidy = 2, genclone = FALSE, sep = ";")
comb_SB1
# ==> Use the same code for SB2, this time sampling 12 individuals without replacement.
```

### Bias identification : Comparison of the 3 methods

In order to have enough samples and because the structure analysis
showed no strong structure between sites SB1 and SB2, we’ll combine the
two sites in one.

Datasets used for the analysis :

``` r
# Zygotes data (n=24)
zyg<-read.genepop(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                           "zygotes_SB_1pop.gen"), ncode = 3L, quiet = TRUE)
zyg

# Duplicated data (n=32)
dup <- read.genepop(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                         "duplicated_SB_1pop.gen"), ncode = 2L, quiet = TRUE)
dup

# Combined (n=15)
comb <- read.genepop(here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", 
                         "combined_SB_1pop.gen"), ncode = 2L, quiet = TRUE)
comb
```

The statistical plan is the following :

- 1 observation = 1 Ne estimation (+ Jackknife CI)
- 1 variable = Ne
- 3 factors : zyg, dup, comb
- 10 replicates for each factor –\> subsample each group 10 times
  (sample 10 individuals without replacement), except for the
  “duplicated data”, were we sample 20 individuals.

``` r
# Individuals sampling

output.dir <- here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", "Replicates_duplicated_new", "Tub")

# Define the number of replicates
n_rep <- 10

# Initialize a list to store the resulting genind objects
sampled_geninds <- vector("list", length = n_rep)

# Perform the sampling for each replicate
for (i in seq_len(n_rep)) {
  
  # Determine the total number of individuals in the genind object
  n_ind <- nrow(dup@tab)
  
  # Set the seed for reproducibility within this replicate
  set.seed(i * 123)
  
  # Sample 50 random individuals without replacement
  sampled_indices <- sample(n_ind, size = 20, replace = FALSE)
  
  # Subset the genind object based on the sampled indices
  sampled_geninds[[i]] <- dup[sampled_indices, ]
}

for (i in seq_along(sampled_geninds)) {
  tmp_file <- tempfile(pattern = "genepop.", fileext = ".txt")
  genind_to_genepop(sampled_geninds[[i]], output = tmp_file)
  filename <- paste0(output.dir, "sampled_20dup_", i, ".txt")
  file.copy(tmp_file, filename, overwrite = TRUE)
  unlink(tmp_file)
}

# The same code is used for the 3 methods, just replacing "duplicated" or "dup" by "zygotes/zyg" or "combined/comb"
```

#### Genetic metrics : Allelic richness, He and Fis

For these 3 methods, the allelic richness, the expected heterozygosity
(He) and the inbreeding coefficient (Fis) are calculated across the 10
replicates.

``` r
# get back the 10 replicates of each group
input_dir <- here("Data", "Data_processed", "Tuber_melanosporum_Taschen2016", "comparison", "Replicates_zygotes")

# Get a list of all Genepop files in the directory
genepop_files <- list.files(path = input_dir, pattern = "*.gen$", full.names = TRUE)

# Read each Genepop file into a genind object and store them as a list
genind_list <- lapply(genepop_files, function(x) {
  genind_obj <- read.genepop(file = x, ncode = 3L, quiet = TRUE)
  attr(genind_obj, "comments") <- ""
  genind_obj
})

zyg <- genind_list

# Do the same for each group

# Calculate He and Fis for each sample - using the package "hierfstat"

# Zygotes
zyg_stats <- lapply(zyg, basic.stats); zyg_stats

# Extract 'Fis' values from '$overall' lists
fis_values <- sapply(zyg_stats, \(x) x[["overall"]][["Fis"]])

# Create a data frame
zyg_Fis <- data.frame("Grp"=c("zyg"), "Fis" = fis_values)

zyg_summ <- lapply(zyg, summary); zyg_summ

# Extract 'Hexp' values from each list
he_values <- t(data.frame(sapply(zyg_summ, \(x) x[["Hexp"]])))

# Create a data frame
zyg_He <- data.frame("Grp"=c("zyg"), "He" = rowMeans(he_values))

# Duplicated
dup_stats <- lapply(dup_new, basic.stats); dup_stats

# Extract 'Fis' values from '$overall' lists
fis_values <- sapply(dup_stats, \(x) x[["overall"]][["Fis"]])

# Create a data frame
dup_Fis <- data.frame("Grp"=c("dup"), "Fis" = fis_values)

dup_summ <- lapply(dup_new, summary); dup_summ

# Extract 'Hexp' values from each list
he_values <- t(data.frame(sapply(dup_summ, \(x) x[["Hexp"]])))

# Create a data frame
dup_He <- data.frame("Grp"=c("dup"),"He" = rowMeans(he_values))


# Combined
comb_stats <- lapply(comb, basic.stats); comb_stats

# Extract 'Fis' values from '$overall' lists
fis_values <- sapply(comb_stats, \(x) x[["overall"]][["Fis"]])

# Create a data frame
comb_Fis <- data.frame("Grp"=c("comb"),"Fis" = fis_values)

comb_summ <- lapply(comb, summary); comb_summ

# Extract 'Hexp' values from each list
he_values <- t(data.frame(sapply(comb_summ, \(x) x[["Hexp"]])))

# Create a data frame
comb_He <- data.frame("Grp"=c("comb"), "He" = rowMeans(he_values))

# Combine all data in a df
He <- dplyr::bind_rows(zyg_He, dup_He, comb_He)
Fis <- dplyr::bind_rows(zyg_Fis, dup_Fis, comb_Fis)

metrics <- dplyr::bind_cols(He, Fis)
metrics <- metrics[,c(1,2,4)]

# Export table to csv
write_csv(metrics, here("Outputs", "Brief_results", "lifecycle_Tub_mel", "Report_analysis", "new_metrics_He_Fis.csv"))
```

#### Ne estimation

The estimation of Ne is performed using NeEstimator, on the 10
replicates of each group (zygotes, duplicated, combined).

## Effect of clonality on Ne estimation

Dataset used : *Boletus edulis* from Hoffman et al. 2020

In this part, we are looking for an effect of the presence of clones in
a genetic dataset on the estimation of Ne.

### Filtering step

``` r
# Original dataset
data <- read.genepop(here("Data", "Data_processed", "Boletus_edulis_Hoffman2020","Boletus_edulis.gen"), ncode = 2L, quiet = TRUE) 

# Missing data correction
# Remove loci with more than 20% of missing data
missing.loc<-missingno(data, type = "loci", cutoff = 0.20);missing.loc

# Remove individuals with more than 20% of missing data
missing.ind<-missingno(missing.loc, type = "genotype", cutoff = 0.20);missing.ind

# Make this cleaned dataset the new main dataset
data<-missing.ind;data

# Correction : Make sample indentifiers consecutive
# Extract unique and ordered rownames (identifier labels)
ordered_ids <- sort(unique(rownames(data@tab)))
# Create a mapping data frame linking old and new identifiers
mapping_df <- data.frame(old_id = ordered_ids, new_id = sequence(length(ordered_ids)), stringsAsFactors = FALSE)
# Apply the mapping to update the rownames attribute
rownames(data@tab) <- mapping_df$new_id

# Structure analysis to find homogenous pops
# Find the best number of groups
grp<-find.clusters(data, max.n.clust = 15) # max.n.clust can be the number of sites or "populations"
# Specify the values according to the plots
### number of PC --> all
### number of clusters --> ideally, the lowest BIC, in practice, the "elbow"

# DAPC 
dapc1 <- dapc(data, pop=grp$grp) # dapc using the groups defined with `find.clusters`
# Specify the values according to the plots

# Display percentages of variation explained by each axis
pourcentages_axes <- dapc1$eig/sum(dapc1$eig) * 100
pourcentages_axes
barplot(pourcentages_axes, main = "Percentage of variation explained by each axis of DAPC",
        xlab = "Axes", ylab = "Percentage of variation")

# Plot the result of DAPC
scatter(dapc1)
scatter(dapc1, posi.da="bottomright",  posi.pca="bottomleft", 
        scree.da = TRUE, bg="white", pch = 20, cstar = 0, cex=1.3, clab=0, 
        legend = T, scree.pca=TRUE)

# Assign the individuals to the new groups (clusters)
groups <- as.factor(grp$grp)
levels(groups) <- c("Pop1", "Pop2", "Pop3")
data@pop <- groups

# Separate the genind object in 3 pops
Pops <- seppop(data)
Pop1 <- Pops$Pop1 ; Pop2 <- Pops$Pop2 ; Pop3 <- Pops$Pop3

# Before going further, let's do a quick analysis of MLG, to see which population will have 
# enough "individuals"
poppr(Pop1)
poppr(Pop2)
poppr(Pop3)
# We take the largest population -> Pop1, with N = 60 and MLG = 31
```

### Bias identification : comparison of 5 treatments

The goal is to compare the estimations of Ne for 5 types of populations
:

- complete (no correction),
- without big clones (\>= 7),
- without big and medium clones (\>= 3),
- without any clones (MLG),
- without any similar MLL.

``` r
# Using Rclone to visualise clones in our Population
# Plot
P.tab <- mlg.table(Pop1)

# Convert the genind object into a dataframe that can be used by Rclone
data.df <- genind2df(Pop1);View(data.df)
popvec <- data.df[,1] # pop info (if present)
dataset <- data.df[,2:ncol(data.df)]# dataset
dataset <- as.data.frame(apply(dataset, 2, function(x) replace(x, is.na(x), "0000"))) # replace NAs by 0000
data2 <- convert_GC(dataset, 2) # change number of digits per allele
row.names(data2) <- seq_len(nrow(data2)) # make rownames consecutive integers starting from 1 to avoid bugs
head(data2)

# MLG list
MLGlist <- MLG_list(data2);MLGlist

# Without big clones
# List of samples to remove = all samples from MLGs with >= 7 clones, except one sample
big_clones <- c(3,4,5,6,7,8,9,35,36,37,38,39,40,54,55,56,57,58,59)
all<-c(1:60)
common_indices <- all %in% big_clones # Find common indices between big_clones and all
new_ind <- all[-which(common_indices)] # Remove common indices from all
new_ind

# Create a new genind object without the big clonal individuals
no_large_clone <- Pop1[new_ind, ]
# Check if the genind object is correct
no_large_clone


# Without medium and large clones
# List of samples to remove = all samples from MLGs with >= 3 clones, except one sample
MLGlist
med_clones <- c(3,4,5,6,7,8,9,20,21,23,22,24,29,35,36,37,38,39,40,48,49,54,55,56,57,58,59)
all<-c(1:60)
common_indices <- all %in% med_clones # Find common indices between big_clones and all
new_ind <- all[-which(common_indices)] # Remove common indices from all
new_ind

# Create a new genind object without the medium and large clonal individuals
no_med_clone <- Pop1[new_ind, ]
# Check if the genind object is correct
no_med_clone

# Without any clone
# Extract one individual from each MLG (= remove all clones)
new_ind <- sapply(MLGlist, `[[`, 1)
new_ind <- as.numeric(new_ind);new_ind
# Create a new genind object without the clonal individuals
noclones <- Pop1[new_ind, ]
# Check if the genind object is correct
noclones

# Without any similar MLL
# MLL analysis

#genetic distances computation, distance on allele differences:
respop <- genet_dist(data2);respop # max 10 differences between genotypes
ressim <- genet_dist_sim(data2, nbrepeat = 1000) #theoretical distribution : sexual reproduction
ressimWS <- genet_dist_sim(data2, genet = TRUE, nbrepeat = 1000) #idem, without selfing

#graph prep.:
p1 <- hist(respop$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), main = "popsim", 
           xlab = "Genetic distances", breaks = seq(0, max(respop$distance_matrix)+1, 1))

p2 <- hist(ressim$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), main = "popSR",
           xlab = "Genetic distances", breaks = seq(0, max(ressim$distance_matrix)+1, 1))

p3 <- hist(ressimWS$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3),
           main = "popSRWS", xlab = "Genetic distances", breaks = seq(0, max(ressimWS$distance_matrix)+1, 1))

limx <- max(max(respop$distance_matrix), max(ressim$distance_matrix), max(ressimWS$distance_matrix))
limy <- 0.30

#graph superposition:
plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), ylim = c(0,limy), main = "", xlab = "Genetic distances")
plot(p2, col = rgb(0.7,0.9,1,0.5), freq = FALSE, add = TRUE)
plot(p3, col = rgb(0.9,0.5,1,0.3), freq = FALSE, add = TRUE)

#adding a legend:
leg.txt <- c("original data","simulated data", "without selfing")
col <- c(rgb(0,0.4,1,1), rgb(0.7,0.9,1,0.5), rgb(0.9,0.5,1,0.3))
legend("top", fill = col, leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white")

#determining alpha2 (= highest genetic distance (i.e. max number of differences) between genotypes that we condider for MLL)
# Can be determined graphically (cutoff) or by looking directly at the data :
table(respop$distance_matrix)
# alpha2 = 1

MLLlist <- MLL_generator(data2, alpha2 = 1) # Change alpha2 value
#or (if problem)
# res <- genet_dist(data2, alpha2 = 2)
# MLLlist <- MLL_generator2(potential_clones = res$potential_clones, res_mlg = MLGlist)

MLLlist

# Extract one individual from each MLL 
new_ind <- sapply(MLLlist, `[[`, 1)
new_ind <- as.numeric(new_ind);new_ind

# Create a new genind object without the clonal individuals
no_sim_MLL <- Pop1[new_ind, ]
# Check if the genind object is correct
no_sim_MLL
```

Now that we have these 5 datasets, with respectively 60, 41, 33, 31 and
27 individuals, we need to standardize the samples for comparison. A
random sampling with n = 25 is realised in each dataset, with 10
replicates.

``` r
# Individuals sampling

output.dir <- here("Data", "Data_processed", "Boletus_edulis_Hoffman2020", "n25_samples", "Replicates_noMLL", "Bol")

# Define the number of replicates
n_rep <- 10

# Initialize a list to store the resulting genind objects
sampled_geninds <- vector("list", length = n_rep)

# Perform the sampling for each replicate
for (i in seq_len(n_rep)) {
  
  # Determine the total number of individuals in the genind object
  n_ind <- nrow(no_sim_MLL@tab)
  
  # Set the seed for reproducibility within this replicate
  set.seed(i * 123)
  
  # Sample 50 random individuals without replacement
  sampled_indices <- sample(n_ind, size = 25, replace = FALSE)
  
  # Subset the genind object based on the sampled indices
  sampled_geninds[[i]] <- no_sim_MLL[sampled_indices, ]
}

for (i in seq_along(sampled_geninds)) {
  tmp_file <- tempfile(pattern = "genepop.", fileext = ".txt")
  genind_to_genepop(sampled_geninds[[i]], output = tmp_file)
  filename <- paste0(output.dir, "sampled_noMLL", i, ".txt")
  file.copy(tmp_file, filename, overwrite = TRUE)
  unlink(tmp_file)
}

# Same code for each dataset, just change the name of the output file and genind object
```

Then Ne is estimated in NeEstimator for the 10 replicates of each
treatment. The results are analysed with statistical tests to assess
significant differences between treatments.

``` r
# Compare Ne estimations

# Data (table containing Ne results for 10 replicates * 5 treatments)
data <- read.csv(here("Outputs", "Brief_results", "clonality", "Clone_correction_bol_edu.csv"),
                 header = TRUE, sep = ";", dec =",")

# Subset data per category

nocorr <- filter(data, Corr == "nocorr")
nolarge <- filter(data, Corr == "nolarge")
nomed <- filter(data, Corr == "nomed")
nomlg <- filter(data, Corr == "nomlg")
nomll <- filter(data, Corr == "nomll")

# Data exploration
boxplot(Ne.~Corr, data = data, main = "Ne estimation as a function of the proportion of clones removed",
        xlab = "Clone correction", ylab = "^Ne")
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
# Statistical tests

# 1-way ANOVA

#Type of plan 
tapply(data$Ne.,data$Corr,length)
```

    ##  nocorr nolarge   nomed   nomlg   nomll 
    ##      10      10      10      10      10

``` r
# Balanced plan with 1 fixed factor

# Application conditions verification
AN<-lm(Ne.~Corr,data)
x11();par(mfrow=c(2,2)); plot(AN) #graphically
# or with tests
shapiro.test(AN$res)# Normality
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  AN$res
    ## W = 0.95983, p-value = 0.08748

``` r
leveneTest(AN)# Levene Test
```

    ## Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    ## factor.

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  4  1.5533 0.2032
    ##       45

``` r
dwtest(AN) # Independence
```

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  AN
    ## DW = 2.387, p-value = 0.7935
    ## alternative hypothesis: true autocorrelation is greater than 0

``` r
# Search for influent points
data$residus<-AN$res
grubbs.test(AN$res, type = 10, opposite = FALSE, two.sided = FALSE)# No significant outlier
```

    ## 
    ##  Grubbs test for one outlier
    ## 
    ## data:  AN$res
    ## G.19 = 2.96801, U = 0.81655, p-value = 0.04789
    ## alternative hypothesis: highest value 7.23 is an outlier

``` r
# Application conditions met --> we can do an ANOVA
anova(AN)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Ne.
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Corr       4 600.94 150.235  23.251 1.812e-10 ***
    ## Residuals 45 290.76   6.461                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(AN)
```

    ## 
    ## Call:
    ## lm(formula = Ne. ~ Corr, data = data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -4.590 -1.585 -0.410  1.385  7.230 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   2.6300     0.8038   3.272 0.002056 ** 
    ## Corrnolarge   1.7400     1.1368   1.531 0.132862    
    ## Corrnomed     4.3200     1.1368   3.800 0.000432 ***
    ## Corrnomlg     6.8600     1.1368   6.035 2.77e-07 ***
    ## Corrnomll     9.6600     1.1368   8.498 6.58e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.542 on 45 degrees of freedom
    ## Multiple R-squared:  0.6739, Adjusted R-squared:  0.6449 
    ## F-statistic: 23.25 on 4 and 45 DF,  p-value: 1.812e-10

``` r
# Post hoc test : Tukey (balanced plan)
TUKEY<-TukeyHSD(aov(Ne.~Corr,data));TUKEY
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Ne. ~ Corr, data = data)
    ## 
    ## $Corr
    ##                diff        lwr       upr     p adj
    ## nolarge-nocorr 1.74 -1.4901286  4.970129 0.5485533
    ## nomed-nocorr   4.32  1.0898714  7.550129 0.0037695
    ## nomlg-nocorr   6.86  3.6298714 10.090129 0.0000027
    ## nomll-nocorr   9.66  6.4298714 12.890129 0.0000000
    ## nomed-nolarge  2.58 -0.6501286  5.810129 0.1738193
    ## nomlg-nolarge  5.12  1.8898714  8.350129 0.0004340
    ## nomll-nolarge  7.92  4.6898714 11.150129 0.0000001
    ## nomlg-nomed    2.54 -0.6901286  5.770129 0.1859471
    ## nomll-nomed    5.34  2.1098714  8.570129 0.0002334
    ## nomll-nomlg    2.80 -0.4301286  6.030129 0.1175543

``` r
letters <- multcompLetters(TUKEY$`Corr`[,4]);letters
```

    ## nolarge   nomed   nomlg   nomll  nocorr 
    ##    "ab"    "ac"    "cd"     "d"     "b"

The results were plotted using the ggplot2 package.

``` r
# Define the order of factors
data$condition <- factor(data$Corr, levels = c("nocorr", "nolarge", "nomed", "nomlg", "nomll"))

# Final plot
fig2 <- ggplot(data, aes(x = condition, y = Ne., fill = condition)) +
  geom_boxplot(alpha=0.5)+ 
  scale_fill_brewer(type = "seq", palette = "YlOrRd", direction = -1) +
  labs(title = "",
       x = "Clone correction",
       y = "^Ne",
       fill = "Clone correction") +
  scale_x_discrete(labels = c("nocorr" = "None", "nolarge" = "Large clones corr.", 
                              "nomed" = "Medium/large clones corr.", "nomlg" = "1 ramet/MLG",
                              "nomll" = "1 ramet/MLL")) +
  annotate("text", x = 1, y = 7, label = "a", size=5, fontface="bold") +
  annotate("text", x = 2, y = 14, label = "ab", size=5, fontface="bold") +
  annotate("text", x = 3, y = 15, label = "bc", size=5, fontface="bold") +
  annotate("text", x = 4, y = 16, label = "cd", size=5, fontface="bold") +
  annotate("text", x = 5, y = 19, label = "d", size=5, fontface="bold") +
  annotate("text", x = 1.8, y = 18, label = "1-way ANOVA : p-value *** R2 = 0.6739", size = 3.5) +
  theme_bw() +
  theme(axis.text.x = element_text(face= "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
    guides(x = guide_axis(n.dodge = 2)) +
  guides(fill = "none")
fig2
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Effect of genetic structure on Ne estimation

Dataset used : *Suillus brevipes* from Branco et al. 2017.

In this dataset, filtering steps including missing data management,
clonal correction and structure analysis had already been done by the
authors. The goal of this part is to estimate Ne for different
populations and groups of populations in order to assess a potential
effect of genetic structure. As this dataset concerns SNP data, the
whole analysis was processed using the Genotoul bioinformatics platform.
The following script details how to standardize the sampling of
individuals (n=12) and creating 10 replicates of 2K SNPs samples.

In the following, abbreviations correspond to the following populations
:

- cal = California
- nocal = Inland (all except California)
- mocal = Mountain California
- cocal = Coastal California
- wycol = Wyoming + Colorado
- minwhi = Minnesota + Canada (Whitecourt)

``` bash
cd /usr/local/bioinfo/src
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
```

Once all datasets are created (in VCF format), they are converted into
the GENEPOP format using the software programme PGDSpider. Ne is
estimated using these new files in NeEstimator.

## Effect of pseudo-replication on Ne estimation

Datasets used : *Suillus brevipes* from Branco et al. 2017 and *Suillus
luteus* from Bazzicalupo et al. 2020.

In these two species genotyped with SNPs, we compare Ne estimations from
datasets containing large VS small datasets, i.e. containing 20,000 SNPs
VS 2,000 SNPs.

### *Suillus brevipes*

The analysis was carried out in *Suillus brevipes* on the 6 populations
identified above (see section ). 12 individuals were sampled in each
population, and 10 subsets were made with 2,000 SNPs and 20,000 SNPs
from the original dataset.

#### 2K SNPs subsets

``` bash
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
```

#### 20K SNPs subsets

``` bash
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
```

#### Ne estimation

Once all datasets are created (in VCF format), they are converted into
the GENEPOP format using the software programme PGDSpider. Ne is
estimated using these new files in NeEstimator.

We want to compare 2 means : mean Ne for 20 000 SNPs and mean Ne for 2
000 SNPs. As we have sampled SNPs randomly from the same pool of SNPs,
the samples are not independent. Then, we’ll use a Wilcoxon test for
paired samples (non parametric).

``` r
# Data
data <- read.csv(here("Outputs", "Brief_results", "pseudorep_Suillus", "Report_analysis", "Suillus_pseudorep_new.csv"), 
                 header = TRUE, sep = ";", dec = ",")

# Subset data per population

all <- filter(data, Pop == "all")
cal <- filter(data, Pop == "cal")
nocal <- filter(data, Pop == "nocal")
mocal <- filter(data, Pop == "mocal")
wycol <- filter(data, Pop == "wycol")
minwhi <- filter(data, Pop == "minwhi")

# Data exploration
boxplot(Ne.~Nbr_SNPs, data = all, main = "Ne estimation from 20K SNPs and 2K SNPs samples",
        xlab = "Number of SNPs", ylab = "^Ne")
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# Wilcoxon test
wilcox.test(Ne.~Nbr_SNPs, data = all, paired=TRUE, alternative="less")
```

    ## 
    ##  Wilcoxon signed rank exact test
    ## 
    ## data:  Ne. by Nbr_SNPs
    ## V = 0, p-value = 0.0009766
    ## alternative hypothesis: true location shift is less than 0

``` r
# Do boxplots and tests for each population by changing data = "..."
```

We can represent the results in a barplot.

``` r
# Define palette
palette_Suillus <- c("#5DA5DA", "#059748")

# Define the order of factors
data$condition <- factor(data$Pop, levels = c("all", "cal", "nocal", "mocal", "minwhi", "wycol"))

# Final plot 2

fig <- ggplot(data, aes(x = condition, y = Ne., fill = Nbr_SNPs)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge", width = 0.7) +
  scale_fill_manual(values = palette_Suillus, labels = c("20K" = "20 000", "2K" = "2 000")) +
  scale_x_discrete(labels = c("all" = "All", "cal" = "California", 
                              "minwhi" = "Min/Can", 
                              "mocal" = "Mount. Cal.", "nocal" = "Inland",
                              "wycol" = "Wyo/Col")) +
  labs(title = "",
       x = "Populations",
       y = "^Ne",
       fill = "Number of SNPs") +
  annotate("text", x = 1, y = 21, label = "***", size=5, fontface="bold") +
  annotate("text", x = 2, y = 45, label = "n.s.", size=4, fontface="bold") +
  annotate("text", x = 3, y = 15, label = "**", size=5, fontface="bold") +
  annotate("text", x = 4, y = 74, label = "n.s.", size=4, fontface="bold") +
  annotate("text", x = 5, y = 41, label = "**", size=5, fontface="bold") +
  annotate("text", x = 6, y = 11, label = "**", size=5, fontface="bold") +
  geom_segment(aes(x = 0.80, xend = 1.20, # all
                   y = 18, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 0.80, xend = 0.80,
                   y = 15, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.20, xend = 1.20,
                   y = 15, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.80, xend = 2.20, # cal
                   y = 41, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.80, xend = 1.80, 
                   y = 38, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.20, xend = 2.20, 
                   y = 38, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.80, xend = 5.20, # minwhi
                   y = 38, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.80, xend = 4.80, 
                   y = 35, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.20, xend = 5.20, 
                   y = 35, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.80, xend = 4.20, # mocal
                   y = 70, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.80, xend = 3.80, 
                   y = 67, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.20, xend = 4.20, 
                   y = 67, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.80, xend = 3.20, # Inland
                   y = 12, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.80, xend = 2.80, 
                   y = 9, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.20, xend = 3.20, 
                   y = 9, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.80, xend = 6.20, # wycol
                   y = 8, yend = 8),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.80, xend = 5.80, 
                   y = 5, yend = 8),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 6.20, xend = 6.20, 
                   y = 5, yend = 8),
               color = "black", linewidth = 0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(r = 50), face= "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        legend.position = "top",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold"))

fig
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### *Suillus luteus*

The analysis was carried out in *Suillus brevipes* on the only
population identified by the authors in their dataset. All individuals
were used (n = 34), and 10 subsets were made with 2,000 SNPs and 20,000
SNPs from the original dataset.

#### 2K SNPs subsets

``` bash
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
```

#### 2K SNPs subsets

``` bash
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
```

#### Ne estimation

Once all datasets are created (in VCF format), they are converted into
the GENEPOP format using the software programme PGDSpider. Ne is
estimated using these new files in NeEstimator.

We want to compare 2 means : mean Ne for 20 000 SNPs and mean Ne for 2
000 SNPs. As we have sampled SNPs randomly from the same pool of SNPs,
the samples are not independent. Then, we’ll use a Wilcoxon test for
paired samples (non parametric).

``` r
# Data
data <- read.csv(here("Outputs", "Brief_results", "pseudorep_Suillus", "Report_analysis", "Pseudorep_Sui_lut.csv"), 
                 header = TRUE, sep = ";", dec = ",")

# There is only 1 population

# Data exploration
boxplot(Ne.~Nbr_SNPs, data = data, main = "Ne estimation from 20K SNPs and 2K SNPs samples",
        xlab = "Number of SNPs", ylab = "^Ne")
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# Wilcoxon test
wilcox.test(Ne.~Nbr_SNPs, data = data, paired=TRUE, alternative="less")
```

    ## 
    ##  Wilcoxon signed rank exact test
    ## 
    ## data:  Ne. by Nbr_SNPs
    ## V = 0, p-value = 0.0009766
    ## alternative hypothesis: true location shift is less than 0

We can represent the results using a boxplot, but it’s not really
relevant.

``` r
# Define palette
palette_Suillus <- c("#5DA5DA", "#059748")

# Plot 

fig_test <- ggplot(data, aes(x = Nbr_SNPs, y = Ne., fill = Nbr_SNPs)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = palette_Suillus) +
  labs(title = "",
       x = "Number of SNPs",
       y = "^Ne") +
  annotate("text", x = 1.5, y = 6750, label = "***", size=5, fontface="bold") +
  geom_segment(aes(x = 1.4, xend = 1.6,
                   y = 6500, yend = 6500),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.4, xend = 1.4,
                   y = 6400, yend = 6500),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.60, xend = 1.6,
                   y = 6400, yend = 6500),
               color = "black", linewidth = 0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(r = 50), face= "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
  guides(fill = "none")

fig_test
```

### Combine the 2 species in 1 figure

The two species are presented side by side in the report, here is the
code to generate Figure 4.

``` r
data <- read.csv(here("Outputs", "Brief_results", "pseudorep_Suillus", "Report_analysis", "Two_Suillus.csv"), 
                 header = TRUE, sep = ";", dec = ",")

# The two plots will be created separately, then combined using grid.arrange()

brevipes <- filter(data, Species == "brevipes")
luteus <- filter(data, Species == "luteus")

# Define palette
palette_Suillus <- c("#5DA5DA", "#059748")

# Plot 1 : Suillus brevipes
brevipes$condition <- factor(brevipes$Pop, levels = c("all", "cal", "nocal", "mocal", "minwhi", "wycol"))

fig1 <- ggplot(brevipes, aes(x = condition, y = Ne., fill = Nbr_SNPs)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge", width = 0.7) +
  scale_fill_manual(values = palette_Suillus, labels = c("20K" = "20 000", "2K" = "2 000")) +
  scale_x_discrete(labels = c("all" = "All", "cal" = "California", 
                              "minwhi" = "Min/Can", 
                              "mocal" = "Mount. Cal.", "nocal" = "Inland",
                              "wycol" = "Wyo/Col")) +
  labs(title = "",
       x = "Suillus brevipes",
       y = "^Ne",
       fill = "Number of SNPs") +
  annotate("text", x = 1, y = 21, label = "***", size=5, fontface="bold") +
  annotate("text", x = 2, y = 45, label = "n.s.", size=4, fontface="bold") +
  annotate("text", x = 3, y = 15, label = "**", size=5, fontface="bold") +
  annotate("text", x = 4, y = 74, label = "n.s.", size=4, fontface="bold") +
  annotate("text", x = 5, y = 41, label = "**", size=5, fontface="bold") +
  annotate("text", x = 6, y = 11, label = "**", size=5, fontface="bold") +
  geom_segment(aes(x = 0.80, xend = 1.20, # all
                   y = 18, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 0.80, xend = 0.80,
                   y = 17, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.20, xend = 1.20,
                   y = 17, yend = 18),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.80, xend = 2.20, # cal
                   y = 41, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.80, xend = 1.80, 
                   y = 40, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.20, xend = 2.20, 
                   y = 40, yend = 41),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.80, xend = 5.20, # minwhi
                   y = 38, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.80, xend = 4.80, 
                   y = 37, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.20, xend = 5.20, 
                   y = 37, yend = 38),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.80, xend = 4.20, # mocal
                   y = 70, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.80, xend = 3.80, 
                   y = 69, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 4.20, xend = 4.20, 
                   y = 69, yend = 70),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.80, xend = 3.20, # Inland
                   y = 12, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 2.80, xend = 2.80, 
                   y = 11, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 3.20, xend = 3.20, 
                   y = 11, yend = 12),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.80, xend = 6.20, # wycol
                   y = 8, yend = 8),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 5.80, xend = 5.80, 
                   y = 7, yend = 8),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 6.20, xend = 6.20, 
                   y = 7, yend = 8),
               color = "black", linewidth = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(margin = margin(r = 50), face= "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold.italic", size = 14),
        legend.position = "top",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold"))

fig1
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
# Plot 2 : Suillus luteus

fig2 <- ggplot(luteus, aes(x = Pop, y = Ne., fill = Nbr_SNPs)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge", width = 0.7) +
  scale_fill_manual(values = palette_Suillus, labels = c("20K" = "20 000", "2K" = "2 000")) +
  scale_x_discrete(labels = c("belgium" = "Belgium")) +
  labs(title = "",
       x = "Suillus luteus",
       y = "^Ne",
       fill = "Number of SNPs") +
  annotate("text", x = 1, y = 4500, label = "***", size=5, fontface="bold") +
  geom_segment(aes(x = 0.80, xend = 1.20, # all
                   y = 4300, yend = 4300),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 0.80, xend = 0.80,
                   y = 4200, yend = 4300),
               color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.20, xend = 1.20,
                   y = 4200, yend = 4300),
               color = "black", linewidth = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(margin = margin(r = 50), face= "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold.italic", size = 14),
        legend.position = "top",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold")) +
  guides(fill = "none")

fig2
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

``` r
Figure_4 <- grid.arrange(fig1, fig2, nrow = 1, ncol = 2, widths = 2:1)
```

![](Complete_workflow_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->
