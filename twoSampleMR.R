# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
rm(list = ls())
#> source: https://mrcieu.github.io/TwoSampleMR/articles/exposure.html
library(TwoSampleMR)

### Exposure data ####
# exposure data: Reading in from a file
#> The file must have a header with column names corresponding to the columns described above
bmi_file <- system.file("extdata", "bmi.txt", package="TwoSampleMR")
bmi_exp_dat <- read_exposure_data(bmi_file)
#> If the Phenotype column is not provided (as is the case in this example) 
#> then it will assume that the phenotype’s name is simply “exposure”. 
#> This is entered in the exposure column. It can be renamed manually:
bmi_exp_dat$exposure <- "BMI"

#> A number of sources of instruments have already been curated and are available for use.
# remotes::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
#> The NHGRI-EBI GWAS catalog contains a catalog of significant associations obtained from GWASs
# To use the NHGRI-EBI GWAS catalog
data(gwas_catalog)
#> For example, to obtain instruments for body mass index using the Speliotes et al 2010 study
bmi_gwas <-subset(gwas_catalog, 
                  grepl("Speliotes", Author) & Phenotype=="Body mass index")
bmi_exp_dat <- format_data(bmi_gwas, type="exposure")

# The IEU GWAS database contains the entire summary statistics for thousands of GWASs
# To obtain a list and details about the available GWASs
ao <- available_outcomes()

# LD Clumping
#> For standard two sample MR it is important to ensure that the instruments for the exposure are independent
#> Once instruments have been identified for an exposure variable, the IEU GWAS database can be used to perform clumping
#> You can provide a list of SNP IDs, the SNPs will be extracted from 1000 genomes data, LD calculated between them, and amongst those SNPs that have LD R-square above the specified threshold only the SNP with the lowest P-value will be retained
#> Note that for the instruments in the R/MRInstruments package the SNPs are already LD clumped
bmi_exp_dat <- clump_data(bmi_exp_dat)


### Outcome data ####
#> Once instruments for the exposure trait have been specified, those variants need to be extracted from the outcome trait
#> If we want to perform MR of BMI against coronary heart disease, we need to identify the SNPs that influence the BMI, and then extract those SNPs from a GWAS on coronary heart disease
#> Get the Locke et al 2014 instruments for BMI as an example
bmi_exp_dat <- extract_instruments(outcomes='ieu-a-2')
#> We now need to find a suitable GWAS for coronary heart disease.
ao[grepl("heart disease", ao$trait), ]
#> The most recent CARDIOGRAM GWAS is ID number ieu-a-7. We can extract the BMI SNPs from this GWAS
chd_out_dat <- extract_outcome_data(snps=bmi_exp_dat$SNP, outcomes='ieu-a-7')

# LD proxies
#> By default if a particular requested SNP is not present in the outcome GWAS then a SNP (proxy) that is in LD with the requested SNP (target) will be searched for instead
# proxies = TRUE or FALSE (TRUE by default)


### Harmonise data ####
#> This means that the effect of a SNP on the exposure and the effect of that SNP on the outcome must each correspond to the same allele
dat <- harmonise_data(exposure_dat=bmi_exp_dat, outcome_dat=chd_out_dat)
# Dealing with strand issue
#> Recent GWASs typically present the effects of a SNP in reference to the allele on the forward strand 
#> But as reference panels are updated the forward strand sometimes changes, and GWASs from a few years ago aren’t guaranteed to be using forward strand conventions
#> There are three options to harmonising the data.
#> 1. Assume all alleles are presented on the forward strand
#> 2. Try to infer the forward strand alleles using allele frequency information
#> 3. Correct the strand for non-palindromic SNPs, but drop all palindromic SNPs
#> By default, the harmonise_data function uses option 2

# Drop duplicate exposure-outcome summary sets
#> After data harmonisation, users may find that their dataset contains duplicate exposure-outcome summary sets
dat <- power_prune(dat, method=1, dist.outcome="binary")


### Perform MR ####
res <- mr(dat)


### Sensitivity analyses ####
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
#> The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis
mr_pleiotropy_test(dat)
# Single SNP analysis
#> To obtain the MR estimates using each of the SNPs singly
res_single <- mr_singlesnp(dat)
# Leave-one-out analysis
#> It is possible to perform a leave-one-out analysis, where the MR is performed again but leaving out each SNP in turn, to identify if a single SNP is driving the association
res_loo <- mr_leaveoneout(dat)


### Plots ####
# Scatter plot
#> Depict the relationship of the SNP effects on the exposure against the SNP effects on the outcome using a scatter plot.
p1 <- mr_scatter_plot(res, dat)
# Forest plot
#> To compare the MR estimates using the different MR methods against the single SNP teststo compare the MR estimates using the different MR methods against the single SNP tests
p2 <- mr_forest_plot(res_single)
# Leave-one-out plot
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
# Funnel plot
#> Asymmetry in a funnel plot is useful for gauging the reliability of a particular MR analysis
#> Funnel plots can be produced using the single SNP results
p4 <- mr_funnel_plot(res_single)


### Report ####
#> A report can be generated that performs all MR analyses, sensitivity analyses, and plots, 
#> and presents them in a single self-contained html web page, word document, or pdf document
mr_report(dat)
