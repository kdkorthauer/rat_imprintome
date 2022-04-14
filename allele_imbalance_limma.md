Assessing allele imbalance in the rat imprintome
================
Keegan Korthauer
25/03/2022

# Setup

## Load libraries

``` r
library(tidyverse)
library(broom)
library(Biobase)
library(ggplot2)
theme_set(theme_bw())
library(limma)
library(ggfortify)
library(edgeR)
library(pheatmap)
library(UpSetR)
```

## Read in and format data

Read in `.txt` file provided by Julien.

``` r
ase.raw <- read.delim("Book2.txt")
```

Print basic info.

``` r
str(ase.raw)
```

    ## 'data.frame':    73059 obs. of  52 variables:
    ##  $ chr                                    : chr  "chr1" "chr1" "chr1" "chr1" ...
    ##  $ start                                  : int  53395 56565 197017 237908 237908 249305 318092 390868 400255 515556 ...
    ##  $ end                                    : int  56311 60411 211255 241172 243880 269871 329799 396476 409676 519170 ...
    ##  $ strand                                 : chr  "-" "-" "-" "-" ...
    ##  $ name                                   : chr  "XR_001835498.1" "XR_589829.2" "XR_589830.2" "XM_017589816.1" ...
    ##  $ ID                                     : chr  "XR_001835498.1" "XR_589829.2" "XR_589830.2" "XM_017589816.1" ...
    ##  $ ExonLength                             : int  1182 1735 2185 1647 1870 2703 2589 1752 2037 1866 ...
    ##  $ Wistar_blood_RNA_PRJEB23955_rep1.9_RPKM: num  0 0 0.0743 0.2956 0.7811 ...
    ##  $ BW_EB_RNA_rep1_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EB_RNA_rep2_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EB_RNA_rep3_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EB_RNA_rep3_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep1_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep2_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep3_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EB_RNA_rep3_WKY_NCrl_F1540_q255_RPM : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EB_RNA_rep1_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EB_RNA_rep1_F334_N_F1540_q255_RPM   : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EB_RNA_rep2_BN_F1540_q255_RPM       : num  0 0 13.1 0 0 ...
    ##  $ BF_EB_RNA_rep2_F334_N_F1540_q255_RPM   : num  0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep1_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep1_F334_N_F1540_q255_RPM   : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep2_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep2_F334_N_F1540_q255_RPM   : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep3_BN_F1540_q255_RPM       : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EB_RNA_rep3_F334_N_F1540_q255_RPM   : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep1_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep1_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep3_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep3_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep4_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BW_EPC_RNA_rep4_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep1_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep1_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep2_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep2_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep3_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ WB_EPC_RNA_rep3_WKY_NCrl_F1540_q255_RPM: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EPC_RNA_rep1_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EPC_RNA_rep1_F334_N_F1540_q255_RPM  : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EPC_RNA_rep2_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ BF_EPC_RNA_rep2_F334_N_F1540_q255_RPM  : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep1_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep1_F334_N_F1540_q255_RPM  : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep2_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep2_F334_N_F1540_q255_RPM  : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep3_BN_F1540_q255_RPM      : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ FB_EPC_RNA_rep3_F334_N_F1540_q255_RPM  : num  0 0 0 0 0 ...

``` r
length(unique(ase.raw$name))
```

    ## [1] 73057

Pull out sample metadata from sample names.

``` r
pd <- data.frame(Sample = colnames(ase.raw)[-c(1:8)])
rownames(pd) <- pd$Sample
attrlist <- strsplit(pd$Sample, "_")
pd <- pd %>%
  mutate(Cross = sapply(attrlist, function(x) x[1]),
         Tissue =  sapply(attrlist, function(x) x[2]),
         Rep =  sapply(attrlist, function(x) x[4]),
         Strain = sapply(attrlist, function(x) x[5]),
         Allele = ifelse(substr(Cross, 1, 1) == substr(Strain, 1, 1),
                         "Maternal", "Paternal"),
         Rat = paste0(Rep, Cross, Tissue),
         Cross_group = ifelse(Cross %in% c("BF", "FB"), "BF/FB", "BW/WB"),
         Cross_direction = ifelse(Cross %in% c("BF", "BW"), "Forward", "Reverse"))

# change rep4 (BW EPC) to rep2 (in that cross rep2 is missing)
pd <- pd %>%
  mutate(Rep = gsub("rep4", "rep2", Rep))

str(pd)
```

    ## 'data.frame':    44 obs. of  9 variables:
    ##  $ Sample         : chr  "BW_EB_RNA_rep1_BN_F1540_q255_RPM" "BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM" "BW_EB_RNA_rep2_BN_F1540_q255_RPM" "BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM" ...
    ##  $ Cross          : chr  "BW" "BW" "BW" "BW" ...
    ##  $ Tissue         : chr  "EB" "EB" "EB" "EB" ...
    ##  $ Rep            : chr  "rep1" "rep1" "rep2" "rep2" ...
    ##  $ Strain         : chr  "BN" "WKY" "BN" "WKY" ...
    ##  $ Allele         : chr  "Maternal" "Paternal" "Maternal" "Paternal" ...
    ##  $ Rat            : chr  "rep1BWEB" "rep1BWEB" "rep2BWEB" "rep2BWEB" ...
    ##  $ Cross_group    : chr  "BW/WB" "BW/WB" "BW/WB" "BW/WB" ...
    ##  $ Cross_direction: chr  "Forward" "Forward" "Forward" "Forward" ...

``` r
table(pd$Cross, pd$Allele, pd$Rep, pd$Tissue)
```

    ## , ,  = rep1,  = EB
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        1        1
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1
    ## 
    ## , ,  = rep2,  = EB
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        1        1
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1
    ## 
    ## , ,  = rep3,  = EB
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        0        0
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1
    ## 
    ## , ,  = rep1,  = EPC
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        1        1
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1
    ## 
    ## , ,  = rep2,  = EPC
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        1        1
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1
    ## 
    ## , ,  = rep3,  = EPC
    ## 
    ##     
    ##      Maternal Paternal
    ##   BF        0        0
    ##   BW        1        1
    ##   FB        1        1
    ##   WB        1        1

There are 3 replicates for each cross in each tissue, except for the BF
cross, for which there are 2 replicates. So in all, there are 11 samples
for each tissue (3 reps in each FB, BW, WB, 2 reps in BF) for a total of
22 samples. Each of these has a maternal and paternal allele count, for
a total of 44 rows in our dataset.

Remove duplicate feature:

``` r
length(ase.raw$name)
```

    ## [1] 73059

``` r
length(unique(ase.raw$name))
```

    ## [1] 73057

``` r
(tab <- which(table(ase.raw$name) == 2))
```

    ## NM_001014271.1 NM_001170534.2 
    ##           7169          14497

``` r
ase.raw <- ase.raw %>%
  filter(!name %in% names(tab))
```

Remove the following features:

1.  features with max RPKMs \< 1  
2.  those that don’t satisfy: RPKM >= 0.5 (total from maternal and
    paternal allele) in at least 2 out of 5 or 6 reps per cross in
    *either* tissue (i.e. required to have 2 reps each for both crosses
    in *either* EB or EPC)
3.  those not on canonical autosomes

``` r
filt1 <- ase.raw %>% 
  select(-c(1:8)) %>%
  apply(., 1, max) >= 1
sum(filt1)
```

    ## [1] 52625

``` r
filt2 <- ase.raw %>% 
  select(-c(1:4, 6:8)) %>%
  pivot_longer(names_to = "Sample",
               values_to = "RPKM",
               cols = -1) %>%
  left_join(pd) %>%
  select(-Sample, -Strain, -Rat) %>%
  group_by(name, Cross, Cross_group, Tissue, Rep) %>%
  summarize(sumRPKM = sum(RPKM)) %>%
  ungroup() %>%
  group_by(name, Tissue, Cross_group) %>%
  summarize(total = n(),
            pass = sum(sumRPKM >= 0.5) >= 2) %>%
  ungroup() %>%
  group_by(name, Tissue) %>%
  summarize(pass = sum(pass) == 2) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(pass = sum(pass) > 0) %>%
  ungroup() %>%
  left_join(ase.raw, ., by = "name") %>%
  pull(pass)
```

    ## Joining, by = "Sample"

    ## `summarise()` has grouped output by 'name', 'Cross', 'Cross_group', 'Tissue'. You can override using the `.groups` argument.

    ## `summarise()` has grouped output by 'name', 'Tissue'. You can override using the `.groups` argument.

    ## `summarise()` has grouped output by 'name'. You can override using the `.groups` argument.

``` r
sum(filt2)
```

    ## [1] 26896

``` r
filt3 <- ase.raw$chr %in% paste0("chr", 1:20)
sum(filt3)
```

    ## [1] 70129

``` r
filt <- filt1 & filt2 & filt3
sum(filt)
```

    ## [1] 26356

``` r
ase <- ase.raw %>% filter(filt)
dim(ase)
```

    ## [1] 26356    52

Let’s put this metadata into a single DGEList object:

``` r
dge <- DGEList(counts = ase %>% select(-c(1:8)), 
               samples = pd,
               genes = ase %>% select(c(1:7)))
dge
```

    ## An object of class "DGEList"
    ## $counts
    ##   BW_EB_RNA_rep1_BN_F1540_q255_RPM BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM
    ## 1                          0.00000                               50.15189
    ## 2                          9.62018                                0.00000
    ## 3                          0.00000                                0.00000
    ## 4                          0.00000                               13.18501
    ## 5                          0.00000                               16.45685
    ##   BW_EB_RNA_rep2_BN_F1540_q255_RPM BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM
    ## 1                                0                             130.226003
    ## 2                                0                               3.990941
    ## 3                                0                               8.459061
    ## 4                                0                              22.774368
    ## 5                                0                              26.591817
    ##   BW_EB_RNA_rep3_BN_F1540_q255_RPM BW_EB_RNA_rep3_WKY_NCrl_F1540_q255_RPM
    ## 1                                0                               68.44488
    ## 2                                0                               13.80127
    ## 3                                0                               18.43286
    ## 4                                0                                0.00000
    ## 5                                0                                0.00000
    ##   WB_EB_RNA_rep1_BN_F1540_q255_RPM WB_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM
    ## 1                         0.000000                              53.918769
    ## 2                         0.000000                               6.549642
    ## 3                         1.058528                              19.185819
    ## 4                         0.000000                               5.424956
    ## 5                         0.000000                              12.371545
    ##   WB_EB_RNA_rep2_BN_F1540_q255_RPM WB_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM
    ## 1                                0                              90.230385
    ## 2                                0                               0.000000
    ## 3                                0                              18.236044
    ## 4                                0                              21.987740
    ## 5                                0                               4.416541
    ##   WB_EB_RNA_rep3_BN_F1540_q255_RPM WB_EB_RNA_rep3_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.0000                              61.156882
    ## 2                           0.0000                               0.000000
    ## 3                          10.5658                               0.000000
    ## 4                           0.0000                               1.467472
    ## 5                           0.0000                               6.273443
    ##   BF_EB_RNA_rep1_BN_F1540_q255_RPM BF_EB_RNA_rep1_F334_N_F1540_q255_RPM
    ## 1                          0.00000                             27.66315
    ## 2                          0.00000                              0.00000
    ## 3                         19.16594                              0.00000
    ## 4                         18.88270                              0.00000
    ## 5                          0.00000                              1.88827
    ##   BF_EB_RNA_rep2_BN_F1540_q255_RPM BF_EB_RNA_rep2_F334_N_F1540_q255_RPM
    ## 1                                0                             62.98712
    ## 2                                0                              0.00000
    ## 3                                0                              0.00000
    ## 4                                0                              0.00000
    ## 5                                0                              0.00000
    ##   FB_EB_RNA_rep1_BN_F1540_q255_RPM FB_EB_RNA_rep1_F334_N_F1540_q255_RPM
    ## 1                                0                            20.427281
    ## 2                                0                             7.647652
    ## 3                                0                             0.000000
    ## 4                                0                             8.251414
    ## 5                                0                            27.471224
    ##   FB_EB_RNA_rep2_BN_F1540_q255_RPM FB_EB_RNA_rep2_F334_N_F1540_q255_RPM
    ## 1                                0                             38.47130
    ## 2                                0                              0.00000
    ## 3                                0                              0.00000
    ## 4                                0                             15.97070
    ## 5                                0                             10.46355
    ##   FB_EB_RNA_rep3_BN_F1540_q255_RPM FB_EB_RNA_rep3_F334_N_F1540_q255_RPM
    ## 1                          0.00000                             17.71398
    ## 2                          0.00000                              0.00000
    ## 3                          0.00000                              0.00000
    ## 4                         11.39601                             10.27409
    ## 5                          0.00000                             11.92741
    ##   BW_EPC_RNA_rep1_BN_F1540_q255_RPM BW_EPC_RNA_rep1_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.00000                                 0.00000
    ## 2                          10.50968                                47.07916
    ## 3                          37.21285                                40.10838
    ## 4                           0.00000                                 0.00000
    ## 5                           0.00000                                 0.00000
    ##   BW_EPC_RNA_rep3_BN_F1540_q255_RPM BW_EPC_RNA_rep3_WKY_NCrl_F1540_q255_RPM
    ## 1                          0.000000                                 0.00000
    ## 2                         17.978514                                21.57422
    ## 3                         38.616862                                23.19968
    ## 4                          0.000000                                 0.00000
    ## 5                          1.477686                                 0.00000
    ##   BW_EPC_RNA_rep4_BN_F1540_q255_RPM BW_EPC_RNA_rep4_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.00000                                 0.00000
    ## 2                          35.61108                                29.32677
    ## 3                          84.06411                                80.55759
    ## 4                           0.00000                                 0.00000
    ## 5                           0.00000                                 0.00000
    ##   WB_EPC_RNA_rep1_BN_F1540_q255_RPM WB_EPC_RNA_rep1_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.00000                                8.771700
    ## 2                           0.00000                                7.346856
    ## 3                          15.76234                               17.187190
    ## 4                           0.00000                                0.000000
    ## 5                           0.00000                                0.000000
    ##   WB_EPC_RNA_rep2_BN_F1540_q255_RPM WB_EPC_RNA_rep2_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.00000                                 0.00000
    ## 2                          17.78305                                 9.17231
    ## 3                          74.50170                                16.70671
    ## 4                           0.00000                                 0.00000
    ## 5                           0.00000                                 0.00000
    ##   WB_EPC_RNA_rep3_BN_F1540_q255_RPM WB_EPC_RNA_rep3_WKY_NCrl_F1540_q255_RPM
    ## 1                           0.00000                                 0.00000
    ## 2                           0.00000                                28.34611
    ## 3                          27.43028                                17.85946
    ## 4                           0.00000                                 0.00000
    ## 5                           0.00000                                 0.00000
    ##   BF_EPC_RNA_rep1_BN_F1540_q255_RPM BF_EPC_RNA_rep1_F334_N_F1540_q255_RPM
    ## 1                          7.273532                               0.00000
    ## 2                          7.758451                              49.13684
    ## 3                         16.405853                              16.24422
    ## 4                          0.000000                               0.00000
    ## 5                          0.000000                               0.00000
    ##   BF_EPC_RNA_rep2_BN_F1540_q255_RPM BF_EPC_RNA_rep2_F334_N_F1540_q255_RPM
    ## 1                           0.00000                               0.00000
    ## 2                          14.00756                              34.88022
    ## 3                          13.86889                              13.93824
    ## 4                           0.00000                               0.00000
    ## 5                           0.00000                               0.00000
    ##   FB_EPC_RNA_rep1_BN_F1540_q255_RPM FB_EPC_RNA_rep1_F334_N_F1540_q255_RPM
    ## 1                                 0                              0.000000
    ## 2                                 0                             44.247491
    ## 3                                 0                             19.900451
    ## 4                                 0                              0.000000
    ## 5                                 0                              1.457906
    ##   FB_EPC_RNA_rep2_BN_F1540_q255_RPM FB_EPC_RNA_rep2_F334_N_F1540_q255_RPM
    ## 1                                 0                              0.000000
    ## 2                                 0                             24.722240
    ## 3                                 0                              9.750985
    ## 4                                 0                              0.000000
    ## 5                                 0                              0.000000
    ##   FB_EPC_RNA_rep3_BN_F1540_q255_RPM FB_EPC_RNA_rep3_F334_N_F1540_q255_RPM
    ## 1                           0.00000                               0.00000
    ## 2                          19.38482                              11.24318
    ## 3                          19.38482                               0.00000
    ## 4                           0.00000                               0.00000
    ## 5                           0.00000                               0.00000
    ## 26351 more rows ...
    ## 
    ## $samples
    ##                                        group lib.size norm.factors
    ## BW_EB_RNA_rep1_BN_F1540_q255_RPM           1 12412628            1
    ## BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM     1 11397366            1
    ## BW_EB_RNA_rep2_BN_F1540_q255_RPM           1 13042735            1
    ## BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM     1 12016808            1
    ## BW_EB_RNA_rep3_BN_F1540_q255_RPM           1 11887085            1
    ##                                                                        Sample
    ## BW_EB_RNA_rep1_BN_F1540_q255_RPM             BW_EB_RNA_rep1_BN_F1540_q255_RPM
    ## BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM
    ## BW_EB_RNA_rep2_BN_F1540_q255_RPM             BW_EB_RNA_rep2_BN_F1540_q255_RPM
    ## BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM
    ## BW_EB_RNA_rep3_BN_F1540_q255_RPM             BW_EB_RNA_rep3_BN_F1540_q255_RPM
    ##                                        Cross Tissue  Rep Strain   Allele
    ## BW_EB_RNA_rep1_BN_F1540_q255_RPM          BW     EB rep1     BN Maternal
    ## BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM    BW     EB rep1    WKY Paternal
    ## BW_EB_RNA_rep2_BN_F1540_q255_RPM          BW     EB rep2     BN Maternal
    ## BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM    BW     EB rep2    WKY Paternal
    ## BW_EB_RNA_rep3_BN_F1540_q255_RPM          BW     EB rep3     BN Maternal
    ##                                             Rat Cross_group Cross_direction
    ## BW_EB_RNA_rep1_BN_F1540_q255_RPM       rep1BWEB       BW/WB         Forward
    ## BW_EB_RNA_rep1_WKY_NCrl_F1540_q255_RPM rep1BWEB       BW/WB         Forward
    ## BW_EB_RNA_rep2_BN_F1540_q255_RPM       rep2BWEB       BW/WB         Forward
    ## BW_EB_RNA_rep2_WKY_NCrl_F1540_q255_RPM rep2BWEB       BW/WB         Forward
    ## BW_EB_RNA_rep3_BN_F1540_q255_RPM       rep3BWEB       BW/WB         Forward
    ## 39 more rows ...
    ## 
    ## $genes
    ##    chr   start     end strand           name             ID ExonLength
    ## 1 chr1 1101665 1120340      - XM_006227600.3 XM_006227600.3       4899
    ## 2 chr1 1181300 1200526      + XM_017590411.1 XM_017590411.1        781
    ## 3 chr1 1207153 1220928      + XR_001835499.1 XR_001835499.1       5555
    ## 4 chr1 1390045 1390532      -    XR_589836.1    XR_589836.1        344
    ## 5 chr1 1390548 1395730      - XM_006227603.3 XM_006227603.3       1130
    ## 26351 more rows ...

# EDA

## Principal components analysis

Let’s plot some PCA plots to see where major variation lies.

``` r
pca_res <- prcomp(t(dge$counts), scale. = TRUE)

autoplot(pca_res, data = dge$samples, colour = 'Tissue')
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
autoplot(pca_res, data = dge$samples, colour = 'Cross')
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

PC 1 is clearly Cross (BF/FB vs BW/WB). PC 2 is clearly tissue.

Let’s look at the next couple PCs.

``` r
autoplot(pca_res, data = dge$samples, colour = 'Cross', x = 3, y = 4)
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
autoplot(pca_res, data = dge$samples, colour = 'Strain', x = 3, y = 4)
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
autoplot(pca_res, data = dge$samples, colour = 'Tissue', x = 3, y = 4)
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

PC 3 looks like strain effects (which strain is the allele from), and PC
4 looks like interaction between cross and tissue.

## Sample-sample correlation heatmap

Here we’ll look at a sample-sample correlation heatmap, which will
further let us look at major sources of variation and spot any potential
outlier samples.

``` r
# Set up color scheme for heatmaps 
bcols<-colorRampPalette(c("#000000" ,"#800000" ,"#FF8000" ,"#FFFF00", "#FFFFFF"))(20)

cc <- data.frame(cor(dge$counts), 
                 row.names = colnames(dge))
range(cc, na.rm=T)
```

    ## [1] 0.2235348 1.0000000

``` r
annot_df <- dge$samples %>%
  select(Cross, Tissue, Allele, Strain)
pheatmap(cc, color = bcols, 
         border_color = NA, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         annotation_col = annot_df, 
         main="Sample-Sample Correlation")
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

A few observations:

-   As expected, tissue is the major source of variation  
-   Cross is the next major source of variation (BF/FB vs BW/WB)  
-   Within a tissue and reciprocal cross, expression is more similar for
    alleles coming from the same strain than from the same type of
    allele (mat/pat)  
-   Within a tissue, reciprocal cross and maternal strain, mat/pat
    alleles are generally more similar to each other than to the
    opposite allele within the same rat

## Allelic expression proportions

Here we want to calculate allelic expression proportions
(Maternal/(Maternal + Paternal)) for each rep.

``` r
ase.long <- pivot_longer(data.frame(dge$counts) %>% 
                           mutate(gene = dge$genes$name),
                         cols = -gene,
                         names_to = "Sample",
                         values_to = "RPKM") %>%
  left_join(pd) %>%
  select(-Sample, -Strain) %>%
  pivot_wider(names_from = Allele,
              values_from = "RPKM") %>%
  mutate(prop = Maternal / (Maternal + Paternal)) 
```

    ## Joining, by = "Sample"

``` r
ase.long
```

    ## # A tibble: 579,832 × 10
    ##    gene   Cross Tissue Rep   Rat   Cross_group Cross_direction Maternal Paternal
    ##    <chr>  <chr> <chr>  <chr> <chr> <chr>       <chr>              <dbl>    <dbl>
    ##  1 XM_00… BW    EB     rep1  rep1… BW/WB       Forward              0       50.2
    ##  2 XM_00… BW    EB     rep2  rep2… BW/WB       Forward              0      130. 
    ##  3 XM_00… BW    EB     rep3  rep3… BW/WB       Forward              0       68.4
    ##  4 XM_00… WB    EB     rep1  rep1… BW/WB       Reverse             53.9      0  
    ##  5 XM_00… WB    EB     rep2  rep2… BW/WB       Reverse             90.2      0  
    ##  6 XM_00… WB    EB     rep3  rep3… BW/WB       Reverse             61.2      0  
    ##  7 XM_00… BF    EB     rep1  rep1… BF/FB       Forward              0       27.7
    ##  8 XM_00… BF    EB     rep2  rep2… BF/FB       Forward              0       63.0
    ##  9 XM_00… FB    EB     rep1  rep1… BF/FB       Reverse             20.4      0  
    ## 10 XM_00… FB    EB     rep2  rep2… BF/FB       Reverse             38.5      0  
    ## # … with 579,822 more rows, and 1 more variable: prop <dbl>

Let’s visualize a histogram of maternal allelic proportions per gene.

``` r
ase.long %>% 
  group_by(gene, Tissue, Cross) %>%
  summarize(meanprop = mean(prop, na.rm = TRUE),
            n = sum(!is.na(prop))) %>%
  filter(n > 1) %>%
  ggplot() +
  geom_histogram(aes(meanprop)) +
  facet_grid(Cross ~ Tissue)
```

    ## `summarise()` has grouped output by 'gene', 'Tissue'. You can override using the `.groups` argument.

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

As we expect, most values are close to 0.5.

### Scatterplots of maternal allele proportions by forward/reverse cross

``` r
ase.long %>% 
  group_by(gene, Tissue, Cross, Cross_group, Cross_direction) %>%
  summarize(meanprop = mean(prop, na.rm = TRUE),
            n = sum(!is.na(prop))) %>%
  filter(n > 1) %>%
  pivot_wider(id_cols = c(gene, Tissue, Cross_group), 
              names_from = Cross_direction, values_from = meanprop) %>%
  ggplot() +
  geom_point(aes(x = Forward, y = Reverse), alpha = 0.2, size = 0.5) +
  facet_grid(Cross_group ~ Tissue)
```

    ## `summarise()` has grouped output by 'gene', 'Tissue', 'Cross', 'Cross_group'. You can override using the `.groups` argument.

    ## Warning: Removed 12452 rows containing missing values (geom_point).

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# Limma voom

## Test within tissue

Here we’ll try out using limma voom to estimate allele-specific effects
within each tissue (EB and EPC).

We’ll build one model to represent the entire dataset. Our model here
includes blocking effects for each rat, so that we can estimate maternal
vs paternal allele effects. We’ll include terms to do so separately in
each allele. Note that there’s no need to include separate terms for
cross, due to the presence of the blocking term for rat (and since we
are not interested in looking at differential allele-specific expression
between the crosses).

``` r
mm <- model.matrix(~ 0 + Rat + Tissue:Allele,
                   data = dge$samples)

# drop redundant terms to get to full rank
mm <- mm[,!grepl("AlleleMaternal", colnames(mm))]
colnames(mm)
```

    ##  [1] "Ratrep1BFEB"              "Ratrep1BFEPC"            
    ##  [3] "Ratrep1BWEB"              "Ratrep1BWEPC"            
    ##  [5] "Ratrep1FBEB"              "Ratrep1FBEPC"            
    ##  [7] "Ratrep1WBEB"              "Ratrep1WBEPC"            
    ##  [9] "Ratrep2BFEB"              "Ratrep2BFEPC"            
    ## [11] "Ratrep2BWEB"              "Ratrep2FBEB"             
    ## [13] "Ratrep2FBEPC"             "Ratrep2WBEB"             
    ## [15] "Ratrep2WBEPC"             "Ratrep3BWEB"             
    ## [17] "Ratrep3BWEPC"             "Ratrep3FBEB"             
    ## [19] "Ratrep3FBEPC"             "Ratrep3WBEB"             
    ## [21] "Ratrep3WBEPC"             "Ratrep4BWEPC"            
    ## [23] "TissueEB:AllelePaternal"  "TissueEPC:AllelePaternal"

``` r
# don't normalize for library size (diff is within sample)
y <- voom(dge, mm, plot = TRUE, 
          normalize.method = "none", 
          lib.size = rep(1, ncol(dge)))
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
fit <- lmFit(y, mm)
fit <- eBayes(fit)
```

First we’ll pull out the results for significant allele effects in EB:

``` r
# pull out test for EB
res_EB <- topTable(fit, coef = "TissueEB:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EB)
```

    ##        chr     start       end strand           name             ID ExonLength
    ## 1394  chr1 194953269 194956281      .   MSTRG.1000.1   MSTRG.1000.1       3012
    ## 11733 chr6 133710206 133721191      - XM_017594511.1 XM_017594511.1       9383
    ## 1395  chr1 195032372 195074327      - XR_001835978.1 XR_001835978.1      13741
    ## 11731 chr6 133658936 133691130      +    NR_131064.1    NR_131064.1       1890
    ## 6874  chr4  29795523  29804517      + XM_008762737.1 XM_008762737.1       2403
    ## 6873  chr4  29795452  29804517      +   MSTRG.8729.1   MSTRG.8729.1       2603
    ##           logFC  AveExpr         t      P.Value    adj.P.Val         B
    ## 1394   3.557530 22.15630  17.15999 6.213208e-15 1.637553e-10 20.692383
    ## 11733 -7.963498 22.61760 -13.44914 1.232470e-12 1.624149e-08  9.210961
    ## 1395   2.913677 26.94109  12.33701 7.542000e-12 5.121784e-08 17.225542
    ## 11731 -8.262576 24.26449 -12.31908 7.773234e-12 5.121784e-08  8.598796
    ## 6874   7.372112 25.86902  11.98352 1.376098e-11 7.253689e-08  8.492042
    ## 6873   7.220388 25.92802  11.62063 2.586129e-11 1.136000e-07  8.325136

``` r
sum(res_EB$adj.P.Val < 0.05)
```

    ## [1] 20

There are 20 significant hits in EB

Let’s pull out the results for significant allele effects in EPC:

``` r
# pull out test for EPC
res_EPC <- topTable(fit, coef = "TissueEPC:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EPC)
```

    ##        chr     start       end strand           name             ID ExonLength
    ## 1772  chr1 215744403 215747080      -    NR_027324.1    NR_027324.1       2325
    ## 7545  chr4 119250244 119292643      - XM_006236839.3 XM_006236839.3       5881
    ## 10334 chr5 157282649 157285295      + XM_006239146.2 XM_006239146.2        933
    ## 10335 chr5 157282717 157285294      +    NM_031598.3    NM_031598.3        766
    ## 13960 chr7 145068238 145113501      + XM_008765753.2 XM_008765753.2       4752
    ## 13961 chr7 145068289 145113507      + NM_001108119.1 NM_001108119.1       4623
    ##           logFC  AveExpr         t      P.Value    adj.P.Val        B
    ## 1772  -6.372358 29.31165 -48.10512 2.758686e-25 7.270792e-21 48.04080
    ## 7545  -7.486167 19.79217 -28.03172 8.776912e-20 1.156621e-15 31.55579
    ## 10334 -8.608692 20.22521 -27.14866 1.843259e-19 1.214524e-15 31.59458
    ## 10335 -8.608692 20.22521 -27.14866 1.843259e-19 1.214524e-15 31.59458
    ## 13960 -6.859438 19.63528 -26.41400 3.477562e-19 1.527577e-15 30.25170
    ## 13961 -6.859438 19.63528 -26.41400 3.477562e-19 1.527577e-15 30.25170

``` r
sum(res_EPC$adj.P.Val < 0.05)
```

    ## [1] 950

There are 950 significant hits in EPC.

Look for overlap between the two:

``` r
ix <- which((res_EB %>% filter(adj.P.Val < 0.05) %>% pull(name)) %in% 
            (res_EPC %>% filter(adj.P.Val < 0.05) %>% pull(name)))
length(ix)
```

    ## [1] 12

``` r
res_EB %>%
  filter(adj.P.Val < 0.05) %>%
  slice(ix) %>%
  pull(name)
```

    ##  [1] "XM_017594511.1" "NR_131064.1"    "XM_008762737.1" "MSTRG.8729.1"  
    ##  [5] "MSTRG.8729.2"   "NR_027324.1"    "MSTRG.4373.2"   "MSTRG.1218.1"  
    ##  [9] "MSTRG.10785.1"  "NM_001009617.1" "MSTRG.8805.1"   "NM_001017504.1"

We see that 60% of the significant hits in EB are also significant in
EPC.

## Test within cross and tissue

Here we’ll try out using limma voom to estimate allele-specific effects
within each tissue (EB and EPC) and cross (BF/FB and WB/BW). The
difference between the previous model is that now we’ll look at effects
in each cross group instead of averaged over both.

``` r
mm <- model.matrix(~ 0 + Rat + Tissue:Cross_group:Allele,
                   data = dge$samples)

# drop redundant terms to get to full rank
mm <- mm[,!grepl("AlleleMaternal", colnames(mm))]
colnames(mm)
```

    ##  [1] "Ratrep1BFEB"                              
    ##  [2] "Ratrep1BFEPC"                             
    ##  [3] "Ratrep1BWEB"                              
    ##  [4] "Ratrep1BWEPC"                             
    ##  [5] "Ratrep1FBEB"                              
    ##  [6] "Ratrep1FBEPC"                             
    ##  [7] "Ratrep1WBEB"                              
    ##  [8] "Ratrep1WBEPC"                             
    ##  [9] "Ratrep2BFEB"                              
    ## [10] "Ratrep2BFEPC"                             
    ## [11] "Ratrep2BWEB"                              
    ## [12] "Ratrep2FBEB"                              
    ## [13] "Ratrep2FBEPC"                             
    ## [14] "Ratrep2WBEB"                              
    ## [15] "Ratrep2WBEPC"                             
    ## [16] "Ratrep3BWEB"                              
    ## [17] "Ratrep3BWEPC"                             
    ## [18] "Ratrep3FBEB"                              
    ## [19] "Ratrep3FBEPC"                             
    ## [20] "Ratrep3WBEB"                              
    ## [21] "Ratrep3WBEPC"                             
    ## [22] "Ratrep4BWEPC"                             
    ## [23] "TissueEB:Cross_groupBF/FB:AllelePaternal" 
    ## [24] "TissueEPC:Cross_groupBF/FB:AllelePaternal"
    ## [25] "TissueEB:Cross_groupBW/WB:AllelePaternal" 
    ## [26] "TissueEPC:Cross_groupBW/WB:AllelePaternal"

``` r
# don't normalize for library size (diff is within sample)
y <- voom(dge, mm, plot = TRUE, 
          normalize.method = "none", 
          lib.size = rep(1, ncol(dge)))
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
fit <- lmFit(y, mm)
fit <- eBayes(fit)
```

### EB

First we’ll pull out the results for significant allele effects in EB
for BF/FB:

``` r
# pull out test in EB for BF/FB
res_EB_BF <- topTable(fit, coef = "TissueEB:Cross_groupBF/FB:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EB_BF)
```

    ##       chr     start       end strand           name             ID ExonLength
    ## 1397 chr1 195052255 195097069      -   MSTRG.1001.3   MSTRG.1001.3      23689
    ## 1396 chr1 195052255 195097069      -   MSTRG.1001.2   MSTRG.1001.2      21402
    ## 1395 chr1 195032372 195074327      - XR_001835978.1 XR_001835978.1      13741
    ## 1400 chr1 195074327 195097069      - NM_001270712.1 NM_001270712.1       1844
    ## 1398 chr1 195074327 195096595      -    NM_031117.2    NM_031117.2       1580
    ## 1399 chr1 195074327 195097069      -    NM_130738.2    NM_130738.2       1709
    ##         logFC  AveExpr        t      P.Value    adj.P.Val        B
    ## 1397 4.838540 27.80701 39.69370 8.362726e-22 2.204080e-17 38.79536
    ## 1396 4.836324 27.74668 36.83215 4.167604e-21 5.492068e-17 37.29836
    ## 1395 4.290345 26.94109 26.31870 5.368126e-18 4.716078e-14 30.97210
    ## 1400 8.035983 25.53191 11.68747 7.513302e-11 4.950514e-07 12.93678
    ## 1398 8.243283 25.38895 11.53285 9.659944e-11 5.091950e-07 12.54236
    ## 1399 8.244065 25.45050 11.20998 1.645814e-10 7.229512e-07 12.22226

``` r
sum(res_EB_BF$adj.P.Val < 0.05)
```

    ## [1] 26

There are 26 significant hits in EB for BF/FB.

Next, we’ll pull out the results for significant allele effects in EB
for BW/WB:

``` r
# pull out test in EB for BW/WB
res_EB_BW <- topTable(fit, coef = "TissueEB:Cross_groupBW/WB:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EB_BW)
```

    ##        chr     start       end strand           name             ID ExonLength
    ## 1397  chr1 195052255 195097069      -   MSTRG.1001.3   MSTRG.1001.3      23689
    ## 1396  chr1 195052255 195097069      -   MSTRG.1001.2   MSTRG.1001.2      21402
    ## 1395  chr1 195032372 195074327      - XR_001835978.1 XR_001835978.1      13741
    ## 1394  chr1 194953269 194956281      .   MSTRG.1000.1   MSTRG.1000.1       3012
    ## 11733 chr6 133710206 133721191      - XM_017594511.1 XM_017594511.1       9383
    ## 11731 chr6 133658936 133691130      +    NR_131064.1    NR_131064.1       1890
    ##           logFC  AveExpr         t      P.Value    adj.P.Val         B
    ## 1397   1.687167 27.80701  29.20114 5.930427e-19 1.563023e-14 33.225925
    ## 1396   1.646730 27.74668  27.52845 2.073762e-18 2.732803e-14 32.053079
    ## 1395   2.059172 26.94109  22.07402 2.163047e-16 1.900309e-12 27.380689
    ## 1394   3.654470 22.15630  20.40958 1.105598e-15 7.284785e-12 22.328942
    ## 11733 -8.994401 22.61760 -13.68539 3.592724e-12 1.893797e-08  8.992183
    ## 11731 -9.645413 24.26449 -11.73532 6.954556e-11 3.054904e-07  7.604060

``` r
sum(res_EB_BW$adj.P.Val < 0.05)
```

    ## [1] 24

There are 24 significant hits in EB for BW/WB.

How much overlap between the crosses in EB?

``` r
ix_EB <- which((res_EB_BF %>% filter(adj.P.Val < 0.05) %>% pull(name)) %in% 
               (res_EB_BW %>% filter(adj.P.Val < 0.05) %>% pull(name)))
length(ix_EB)
```

    ## [1] 9

``` r
res_EB_BF %>% 
  filter(adj.P.Val < 0.05) %>%
  slice(ix_EB) %>%
  pull(name)
```

    ## [1] "MSTRG.1001.3"   "MSTRG.1001.2"   "XR_001835978.1" "XM_017594511.1"
    ## [5] "XM_008762737.1" "NR_131064.1"    "MSTRG.8729.1"   "MSTRG.8729.2"  
    ## [9] "NR_027324.1"

We see that 34.6% of the significant EB hits in BF/FB are also
significant in BW/WB.

### EPC

First we’ll pull out the results for significant allele effects in EPC
for BF/FB:

``` r
# pull out test in EPC for BF/FB
res_EPC_BF <- topTable(fit, coef = "TissueEPC:Cross_groupBF/FB:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EPC_BF)
```

    ##        chr     start       end strand           name             ID ExonLength
    ## 1772  chr1 215744403 215747080      -    NR_027324.1    NR_027324.1       2325
    ## 7545  chr4 119250244 119292643      - XM_006236839.3 XM_006236839.3       5881
    ## 13960 chr7 145068238 145113501      + XM_008765753.2 XM_008765753.2       4752
    ## 13961 chr7 145068289 145113507      + NM_001108119.1 NM_001108119.1       4623
    ## 10334 chr5 157282649 157285295      + XM_006239146.2 XM_006239146.2        933
    ## 10335 chr5 157282717 157285294      +    NM_031598.3    NM_031598.3        766
    ##           logFC  AveExpr         t      P.Value    adj.P.Val        B
    ## 1772  -6.078665 29.31165 -35.61024 8.588394e-21 2.263557e-16 37.78614
    ## 7545  -7.822398 19.79217 -20.01444 1.657659e-15 2.184462e-11 21.74387
    ## 13960 -7.196164 19.63528 -19.17291 4.025077e-15 2.652123e-11 20.86244
    ## 13961 -7.196164 19.63528 -19.17291 4.025077e-15 2.652123e-11 20.86244
    ## 10334 -8.789227 20.22521 -18.14240 1.252205e-14 5.500521e-11 20.81660
    ## 10335 -8.789227 20.22521 -18.14240 1.252205e-14 5.500521e-11 20.81660

``` r
sum(res_EPC_BF$adj.P.Val < 0.05)
```

    ## [1] 636

There are 636 significant hits in EPC for BF/FB.

Next, we’ll pull out the results for significant allele effects in EPC
for BW/WB:

``` r
# pull out test in EPC for BW/WB
res_EPC_BW <- topTable(fit, coef = "TissueEPC:Cross_groupBW/WB:AllelePaternal",
                   n = Inf, sort.by = "P")
head(res_EPC_BW)
```

    ##         chr     start       end strand           name             ID ExonLength
    ## 1772   chr1 215744403 215747080      -    NR_027324.1    NR_027324.1       2325
    ## 21615 chr14 100151517 100184192      - NM_001025750.1 NM_001025750.1       1648
    ## 21616 chr14 100157421 100184186      - XM_006251520.1 XM_006251520.1       1256
    ## 5059   chr3  44150493  44225723      - XM_006234202.3 XM_006234202.3       2135
    ## 7545   chr4 119250244 119292643      - XM_006236839.3 XM_006236839.3       5881
    ## 5060   chr3  44150982  44177689      - NM_001012086.1 NM_001012086.1       1673
    ##           logFC  AveExpr         t      P.Value    adj.P.Val        B
    ## 1772  -6.632412 29.31165 -37.64190 2.613981e-21 6.889409e-17 38.94490
    ## 21615 -8.593632 20.01233 -21.64879 3.245755e-16 3.115043e-12 23.80695
    ## 21616 -8.569519 20.00905 -21.55717 3.545731e-16 3.115043e-12 23.74119
    ## 5059  -7.788770 19.92510 -19.41451 3.109095e-15 1.799582e-11 21.98982
    ## 7545  -7.195121 19.79217 -19.18455 3.975050e-15 1.799582e-11 21.46283
    ## 5060  -7.639027 19.89227 -19.15649 4.096787e-15 1.799582e-11 21.72928

``` r
sum(res_EPC_BW$adj.P.Val < 0.05)
```

    ## [1] 613

There are 613 significant hits in EPC for BW/WB.

How much overlap between the crosses in EPC?

``` r
ix_EPC <- which((res_EPC_BF %>% filter(adj.P.Val < 0.05) %>% pull(name)) %in% 
               (res_EPC_BW %>% filter(adj.P.Val < 0.05) %>% pull(name)))
length(ix_EPC)
```

    ## [1] 352

``` r
res_EPC_BF %>% 
  filter(adj.P.Val < 0.05) %>%
  slice(ix_EPC) %>%
  pull(name)
```

    ##   [1] "NR_027324.1"    "XM_006236839.3" "XM_008765753.2" "NM_001108119.1"
    ##   [5] "XM_006239146.2" "NM_031598.3"    "NM_001025138.1" "XM_006234202.3"
    ##   [9] "NM_001012086.1" "XM_006241412.3" "XM_006241411.2" "NM_001130583.1"
    ##  [13] "MSTRG.1216.2"   "XM_008760149.1" "NM_001013194.2" "NM_001310046.1"
    ##  [17] "NM_001310047.1" "NM_001005897.1" "NM_133572.1"    "XM_008769548.1"
    ##  [21] "XM_006235003.1" "NM_001025750.1" "XM_006251520.1" "NM_053455.2"   
    ##  [25] "XM_017598421.1" "XM_017598419.1" "NM_001134602.1" "XM_017598420.1"
    ##  [29] "XM_006236597.3" "XM_001071128.6" "XM_006236596.3" "NM_130421.2"   
    ##  [33] "XM_008758601.1" "NM_001108444.1" "NM_001025638.2" "MSTRG.1856.1"  
    ##  [37] "NM_013096.1"    "NM_001008847.2" "NM_012758.1"    "NM_001127304.1"
    ##  [41] "XM_006237322.3" "XM_017592704.1" "NM_001191896.1" "NM_001003707.1"
    ##  [45] "NM_001037096.1" "XR_001837378.1" "XR_001837377.1" "NM_001276721.1"
    ##  [49] "NM_001033686.1" "XM_017591004.1" "NM_001047918.2" "NM_019336.1"   
    ##  [53] "XM_006249974.3" "NM_012886.2"    "XM_003749638.4" "NM_001008884.2"
    ##  [57] "XM_006241939.3" "XM_008770160.2" "XM_008762737.1" "XM_006241162.3"
    ##  [61] "XR_596699.2"    "MSTRG.8729.2"   "MSTRG.8729.1"   "XR_001841869.1"
    ##  [65] "XM_006254213.3" "XM_006254212.3" "XM_017597720.1" "XM_006254216.3"
    ##  [69] "XM_006254210.3" "XM_006254211.3" "XM_006254215.3" "MSTRG.5575.3"  
    ##  [73] "XM_006254214.3" "XM_001063315.6" "NM_172019.2"    "MSTRG.5575.2"  
    ##  [77] "MSTRG.5575.1"   "MSTRG.5575.4"   "XM_006242140.3" "XR_001837919.1"
    ##  [81] "XR_001837920.1" "NM_001109599.1" "XM_008760160.2" "XM_017589502.1"
    ##  [85] "XM_017589501.1" "NM_001004084.2" "NM_013087.2"    "XM_017600940.1"
    ##  [89] "XM_017594511.1" "NM_133555.1"    "XM_008759821.2" "NM_001008831.3"
    ##  [93] "NM_001170558.1" "NM_001108569.3" "NM_053663.1"    "NM_012523.2"   
    ##  [97] "XM_006230269.3" "NM_001033998.2" "NM_012512.2"    "NM_001170560.1"
    ## [101] "XM_006255028.3" "XM_006255025.3" "XM_017600937.1" "XM_006256061.3"
    ## [105] "XM_006255031.3" "XM_008766341.2" "XM_006255029.3" "XM_017600938.1"
    ## [109] "XM_017600935.1" "XM_017600936.1" "XM_017600939.1" "XM_017595419.1"
    ## [113] "XM_017600934.1" "XM_006255032.3" "XM_017600941.1" "XM_017591207.1"
    ## [117] "NM_016994.2"    "NM_001271215.1" "XM_017591369.1" "XM_017591705.1"
    ## [121] "NM_001107370.1" "NM_022205.3"    "XM_017591206.1" "XM_017591208.1"
    ## [125] "XM_008767630.2" "NM_031140.1"    "NR_131064.1"    "XM_017591755.1"
    ## [129] "XM_017600605.1" "XM_017591704.1" "NM_053634.1"    "XM_006235066.2"
    ## [133] "NM_001107737.1" "XM_008767946.2" "NM_012651.2"    "XM_001059150.6"
    ## [137] "NM_017196.3"    "XM_006236838.3" "NM_130411.2"    "NM_133416.1"   
    ## [141] "XM_017588709.1" "XM_003750140.4" "XM_006240004.3" "XM_017598196.1"
    ## [145] "XM_017600869.1" "XM_006246382.3" "XM_017591706.1" "NM_173153.2"   
    ## [149] "XM_006240005.3" "NM_019337.1"    "NM_001008765.1" "XR_001840564.1"
    ## [153] "XM_006240003.3" "XM_006230827.3" "NM_001108511.1" "XR_001840565.1"
    ## [157] "NM_053313.1"    "MSTRG.5974.1"   "XM_017600872.1" "XM_008772099.2"
    ## [161] "XM_017600868.1" "XM_017600870.1" "XM_017600871.1" "XM_017600866.1"
    ## [165] "XM_017600873.1" "XM_008772100.2" "XR_001840567.1" "NM_001191695.1"
    ## [169] "XM_006230326.3" "MSTRG.5576.1"   "XM_008771159.1" "XM_006252759.2"
    ## [173] "XM_006252762.3" "XM_006230530.2" "NM_001047878.1" "XM_006252760.3"
    ## [177] "XR_001840563.1" "NM_001106314.1" "XM_017600607.1" "XM_017600608.1"
    ## [181] "XM_017600606.1" "XR_001841736.1" "XR_001841737.1" "XM_017600609.1"
    ## [185] "XM_008764215.2" "NM_017226.1"    "NM_138881.1"    "XM_017600867.1"
    ## [189] "XM_006256063.3" "XM_006256062.2" "XM_006236954.3" "XM_006236956.2"
    ## [193] "XM_008763279.2" "NM_053889.1"    "NM_001003691.1" "XM_008769494.2"
    ## [197] "NM_053734.2"    "XM_006230325.3" "NM_001009617.1" "NM_001017504.1"
    ## [201] "XM_017596881.1" "MSTRG.12972.1"  "NM_020542.2"    "XM_017590500.1"
    ## [205] "NM_001012029.1" "XR_593071.1"    "XR_593072.1"    "XM_006242190.3"
    ## [209] "NM_001191789.1" "XM_006234532.3" "XM_006249910.2" "NM_001109890.1"
    ## [213] "NM_138507.2"    "XM_008769534.2" "NM_001109887.1" "NM_001109889.1"
    ## [217] "NM_001109888.1" "XM_006249912.3" "NM_001005892.2" "NM_001129997.1"
    ## [221] "XR_001836904.1" "XM_017596882.1" "MSTRG.2715.1"   "XM_008765689.2"
    ## [225] "XM_017593115.1" "XM_006237924.3" "NM_080479.2"    "NM_001271109.1"
    ## [229] "XR_001841873.1" "NM_001271217.1" "NM_001271110.1" "NM_019169.2"   
    ## [233] "XM_017597823.1" "XR_001841872.1" "XM_008761362.2" "NM_001007602.1"
    ## [237] "XR_001836467.1" "XM_006246770.3" "XM_017597454.1" "NM_017020.3"   
    ## [241] "XM_008761127.2" "NM_001024273.1" "XM_006244178.3" "NM_001106685.1"
    ## [245] "XM_006232518.3" "NM_053372.1"    "NM_001191659.1" "XM_017599193.1"
    ## [249] "NM_012740.3"    "XM_006256065.2" "XM_017591011.1" "NM_001005539.1"
    ## [253] "XM_017598310.1" "XR_001836169.1" "XR_001840566.1" "XM_008768267.2"
    ## [257] "XM_008768266.2" "XM_006250147.3" "NM_170789.2"    "XM_017598309.1"
    ## [261] "NM_001205304.1" "XM_006247499.3" "XM_006247502.3" "XM_006247501.3"
    ## [265] "XM_006247498.3" "XM_006247500.3" "NM_001013062.1" "MSTRG.3132.6"  
    ## [269] "XM_008762963.2" "MSTRG.1218.1"   "XM_017589340.1" "XM_017589339.1"
    ## [273] "NM_001107364.1" "XM_008765663.2" "XR_354304.3"    "XM_008768951.2"
    ## [277] "NM_001164143.3" "XM_008768952.2" "XM_008766055.2" "NM_001107115.1"
    ## [281] "NM_001164142.3" "XM_008767060.2" "XM_017596794.1" "XM_017589341.1"
    ## [285] "XM_017589338.1" "XM_017589337.1" "XM_017589335.1" "XM_017589336.1"
    ## [289] "XM_017589334.1" "XM_017589333.1" "XM_017596795.1" "XM_017590771.1"
    ## [293] "NM_001014843.2" "XM_008762867.2" "XM_006249715.3" "NM_001017381.1"
    ## [297] "NM_031574.1"    "XM_008763369.2" "NM_012843.2"    "NM_153721.1"   
    ## [301] "XM_006256069.3" "XM_017592823.1" "NM_001109247.1" "XM_017590770.1"
    ## [305] "XM_017590772.1" "XM_017590773.1" "XM_017601716.1" "XR_001839132.1"
    ## [309] "XR_001839816.1" "XM_006246453.3" "XM_006230916.2" "XM_006236955.3"
    ## [313] "XM_017597078.1" "XM_002728846.4" "XM_006230915.2" "NM_001106460.1"
    ## [317] "XM_006246455.3" "XM_008762055.2" "XR_001837601.1" "XM_006234627.3"
    ## [321] "MSTRG.8805.1"   "XM_008762058.2" "NM_053519.1"    "XM_008762056.2"
    ## [325] "XR_001839131.1" "XR_001839130.1" "XM_006234629.3" "NM_012588.2"   
    ## [329] "XM_008761107.2" "NM_013171.1"    "XM_017591492.1" "XM_006234631.3"
    ## [333] "XM_017591491.1" "XM_006242277.3" "XM_006240473.3" "XM_017594637.1"
    ## [337] "NM_001077671.1" "NM_001012093.1" "MSTRG.11606.1"  "XM_006242280.3"
    ## [341] "XM_006242278.3" "XM_006242279.3" "XM_006242281.3" "XR_592459.2"   
    ## [345] "NM_001108052.1" "NM_017154.1"    "XM_006240474.3" "XM_017591209.1"
    ## [349] "MSTRG.8582.1"   "NM_001191642.1" "XM_017594813.1" "XM_006242383.3"

We see that 55.3% of the significant EPC hits in BF/FB are also
significant in BW/WB.

### Overlap between crosses and tissues

Any hits across all comparisons (in each cross and tissue)?

``` r
ix_all <- which((res_EB_BF %>% filter(adj.P.Val < 0.05) %>% slice(ix_EB) %>% pull(name)) %in% 
                (res_EPC_BF %>% filter(adj.P.Val < 0.05) %>% slice(ix_EPC) %>% pull(name)))
length(ix_all)
```

    ## [1] 6

``` r
res_EB_BF %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice(ix_EB) %>% 
  slice(ix_all) %>% 
  pull(name)
```

    ## [1] "XM_017594511.1" "XM_008762737.1" "NR_131064.1"    "MSTRG.8729.1"  
    ## [5] "MSTRG.8729.2"   "NR_027324.1"

UpSet Plot to visualize overlaps between four lists of hits.

``` r
listInput <- list(
  EB_BF = res_EB_BF %>% filter(adj.P.Val < 0.05) %>% pull(name),
  EB_BW = res_EB_BW %>% filter(adj.P.Val < 0.05) %>% pull(name),
  EPC_BF = res_EPC_BF %>% filter(adj.P.Val < 0.05) %>% pull(name),
  EPC_BW = res_EPC_BW %>% filter(adj.P.Val < 0.05) %>% pull(name)
)
upset(fromList(listInput), order.by = "freq")
```

![](allele_imbalance_limma_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->
