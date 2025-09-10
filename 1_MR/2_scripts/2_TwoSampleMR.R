#Script to run the Mendelian Randomisation itself via TwoSampleMR and MendelianRandomization
#Author: Jonathan Chan
#Date: 2024-04-30

###--------------------------------------------------------------------------------
library(tidyverse)
library(TwoSampleMR)
library(tictoc)
library(ggmanh)
library(readxl)
library(MRPRESSO)

args = commandArgs(trailingOnly = TRUE) #Allows input of arguments from Rscript

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied", call. = FALSE)
}

exposure_input_path <- args[1]
outcome_input_path <- args[2]
selected_exposure_snps_string <- args[3]
output_path <- args[4]
exposure_name <- args[5]
outcome_name <- args[6]
ld_clumping <- args[7]
ld_eur_bed_file <- args[8]
plink_binary_path <- args[9]

if (ld_clumping == 'TRUE') {
  ld_clumping <- TRUE
} else if (ld_clumping == 'FALSE' | !exists('ld_clumping')) {
  ld_clumping <- FALSE
}

if (selected_exposure_snps_string != 'FALSE') {
  #Assume that the selected exposure SNPs come in as a single string
  selected_exposure_snps <- str_split_1(selected_exposure_snps_string, ',')

  print(selected_exposure_snps)
}

#Rescomp test code for bidirectional MR
# Direction 1
# exposure_input_path <- './levin22_hf_gwas_ssf/final/levin22_hf_gwas_ssf.h.tsv.gz'
# outcome_input_path <- './tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'
# selected_exposure_snps <- '../1_data/gwas_associations/levin22_hf_associations.tsv'
# output_path <- '../3_output/2_MR/'
# exposure_name <- 'HF'
# outcome_name <- 'HCM'
# ld_clumping <- F

#Direction 2
#32 GWS SNPs for HCM for HCM -> HF
# exposure_input_path <- './tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'
# outcome_input_path <- './levin22_hf_gwas_ssf/final/levin22_hf_gwas_ssf.h.tsv.gz'
# selected_exposure_snps <- c("rs1048302","rs11687178","rs3845778","rs2540277","rs6747402","rs7612736","rs4894803","rs2191446","rs11748963","rs3176326","rs12210733","rs66520020","rs7824244","rs35006907","rs2645210","rs2177843","rs11196085","rs17617337","rs12270374","rs182427065","rs7487962","rs41306688","rs113907726","rs8006225","rs8033459","rs28768976","rs7210446","rs2644262","rs6566955","rs12460541","rs62222424","rs5760054")
# output_path <- '../3_output/2_MR/'
# exposure_name <- 'HCM'
# outcome_name <- 'HF'
# ld_clumping <- F

#Checking maxLAVi -> HCM
# exposure_input_path <- './thanaj22_maxLAVi_gwas_ssf/final/thanaj22_maxLAVi_gwas_ssf.h.tsv.gz'
# outcome_input_path <- './tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'
# selected_exposure_snps <- '../1_data/gwas_associations/thanaj22_maxLAVi_associations.tsv'
# output_path <- '../3_output/2_MR/'
# exposure_name <- 'maxLAVi'
# outcome_name <- 'HCM'
# ld_clumping <- F

#---------------------------------------------------------------------------------
#Load the exposure GWAS summary statistics + Extract the instruments
tic()
print('Importing the exposure GWAS summary statistics')

if (any(str_detect(selected_exposure_snps, 'tsv'))) {
  #i.e if a file path is passed instead of a list of SNPs
  selected_exposure_snps <- read_tsv(selected_exposure_snps) %>%
    filter(`P-VALUE` <= 5 * 10^-8) #Filter for only GWS SNPs in the passed in file
  selected_exposure_snps <- selected_exposure_snps$SNPS
} else if (any(str_detect(selected_exposure_snps, 'csv'))) {
  selected_exposure_snps <- read_csv(selected_exposure_snps) %>%
    filter(pValue <= 5 * 10^-8) #Filter for only GWS SNPs in the passed in file
  selected_exposure_snps <- selected_exposure_snps$dbSNP
} else if (any(str_detect(selected_exposure_snps, '.xlsx'))) {
  selected_exposure_snps <- read_excel(selected_exposure_snps, sheet = 'Sheet1')
  selected_exposure_snps <- selected_exposure_snps$rsid
}

#This is the imported summary statisics filtered for those at approx. GWS (i.e pval < 10^-5)
full_exposure_summstats <- read_tsv(exposure_input_path) %>%
  filter(p_value < 10^-5) %>% #Use arbitrary threshold for GWAS significance
  mutate(Phenotype = exposure_name)

#This is as per GWAS-SSF formatting as per EBI GWAS Catalog after harmonisation
exposure_summstats <- format_data(
  full_exposure_summstats,
  type = 'exposure',
  snps = selected_exposure_snps,
  snp_col = 'rsid',
  beta_col = 'beta',
  se_col = 'standard_error',
  eaf_col = 'effect_allele_frequency',
  effect_allele_col = 'effect_allele',
  other_allele_col = 'other_allele',
  pval_col = 'p_value',
  chr_col = 'chromosome',
  pos_col = 'base_pair_location'
)
toc()

#-------------------------------------------------------------------------------
tic()
#LD clumping to only return independent SNPs
if (isTRUE(ld_clumping)) {
  print('LD clumping the associated SNPs with exposure')
  exposure_clumped <- clump_data(
    exposure_summstats,
    bfile = ld_eur_bed_file,
    plink_bin = plink_binary_path
  ) #Uses 1KGP EUR LD Panel

  #Plot a Manhattan to determine validity of instruments after LD clumping
  g2 <- manhattan_plot(
    x = mutate(exposure_clumped, pos.exposure = as.numeric(pos.exposure)),
    pval.colname = "pval.exposure",
    chr.colname = "chr.exposure",
    pos.colname = "pos.exposure",
    plot.title = str_c(
      'Manhattan plot for LD clumped SNPs with pvalue < 10^-5 for phenotype ',
      exposure_name
    ),
    y.label = "-log10(pval)",
    rescale = F,
    label.colname = 'SNP'
  )
  #print(g)
  ggsave(
    str_c(
      output_path,
      exposure_name,
      '_to_',
      outcome_name,
      '/',
      exposure_name,
      '_ldclumped_manhattan.png'
    ),
    g2,
    dpi = 600
  )
  rm(g2)

  if (nrow(exposure_clumped) == 0) {
    print('LD clumping removed all SNPs or you ran out of API calls')
    stop()
  }
} else {
  exposure_clumped <- exposure_summstats
}

toc()

# print(head(exposure_clumped))
exposure_instruments <- exposure_clumped
rm(exposure_clumped, full_exposure_summstats)

if (nrow(exposure_instruments) == 0) {
  print('You have no instruments for the exposure')
  stop()
}

#--------------------------------------------------------------------------------
#Load the outcome GWAS summary statistics
tic()
print('Importing the outcome GWAS summary statistics')

full_outcome_summstats <- read_tsv(outcome_input_path) %>%
  mutate(Phenotype = outcome_name)

outcome_summstats <- format_data(
  full_outcome_summstats,
  snps = exposure_instruments$SNP,
  type = 'outcome',
  snp_col = 'rsid',
  beta_col = 'beta',
  se_col = 'standard_error',
  eaf_col = 'effect_allele_frequency',
  effect_allele_col = 'effect_allele',
  other_allele_col = 'other_allele',
  pval_col = 'p_value',
  chr_col = 'chromosome',
  pos_col = 'base_pair_location'
)

if (nrow(outcome_summstats) != nrow(exposure_instruments)) {
  print(
    'Failed to identify the some of the exposure instruments in the outcome GWAS so ignoring those instruments'
  )
  missing_snps <- exposure_instruments$SNP[
    !exposure_instruments$SNP %in% outcome_summstats$SNP
  ]
  print(str_c(
    'The following ',
    nrow(missing_snps),
    ' SNPs of original ',
    nrow(exposure_instruments),
    ' instruments are not present in the outcome summary statistics: ',
    missing_snps
  ))
} else if (nrow(outcome_summstats) == 0) {
  print(
    'None of the instrument SNPs were found in the outcome summary statistics after formatting'
  )
  stop()
}

toc()
rm(full_outcome_summstats)

#-------------------------------------------------------------------------------
#Harmonise the data
print('Harmonising exposure and outcome data')
harmonised_data <- harmonise_data(
  exposure_instruments,
  outcome_summstats,
  action = 1
) #Assume that all alleles are presented on the forward strand
rm(
  exposure_clumped,
  exposure_instruments,
  exposure_summstats,
  outcome_summstats
)

#--------------------------------------------------------------------------------
#Run MR
print('Running MR')

#Mainline MR analyses
##Select IVW with MRE if n_instruments >5, otherwise fixed-effects as per Burgess et al, 2023
if (nrow(harmonised_data) == 1) {
  mainline_mr_method_list <- c('mr_wald_ratio')
} else {
  mainline_mr_method_list <- c('mr_ivw')
}
mainline_mr_results <- mr(
  harmonised_data,
  method_list = mainline_mr_method_list
) #Runs MR with a bunch of different methods- call mr_method_list() to see\

#Sensitivity analyses
sensitivity_mr_method_list <- c(
  mainline_mr_method_list,
  'mr_weighted_median',
  'mr_weighted_mode',
  'mr_egger_regression'
)
sensitivity_mr_results <- mr(
  harmonised_data,
  method_list = sensitivity_mr_method_list
)
pleiotropy_results <- mr_pleiotropy_test(harmonised_data) #MR Egger test for directional pleiotropy via evaluating MR Egger intercept difference from 0
res_single <- mr_singlesnp(harmonised_data) #Single-SNP analysis
res_loo <- mr_leaveoneout(harmonised_data) #Leave-one-out analysis

#Plotting plots
print('Plotting output plots')
theme_set(theme_classic())
scatter <- mr_scatter_plot(mainline_mr_results, harmonised_data)
scatter2 <- mr_scatter_plot(sensitivity_mr_results, harmonised_data)
# print(scatter[[1]])
forest <- mr_forest_plot(res_single)
#print(forest[[1]])
loo_forest <- mr_leaveoneout_plot(res_loo)
#print(loo_forest[[1]])

#Add an output step to output odds-ratio with 95% confidence interval (instead of log-odds with standard error)
safe_generate_odds_ratios <- safely(generate_odds_ratios)
or_result <- safe_generate_odds_ratios(sensitivity_mr_results)

# Output the harmonised instrument details (Beta and SE for exposure/outcome)
print("Writing individual instrument details (beta/se) to TSV file.")
instrument_details_df <- harmonised_data %>%
  dplyr::select(
    SNP,
    effect_allele.exposure, # Keep alleles for reference
    other_allele.exposure,
    eaf.exposure, # Keep EAF for context
    beta.exposure,
    se.exposure,
    pval.exposure, # Keep p-values
    beta.outcome,
    se.outcome,
    pval.outcome,
    mr_keep # Indicates if SNP was kept after harmonisation checks
  )

instrument_details_filename <- file.path(
  output_path,
  paste0(exposure_name, '_to_', outcome_name, '_instrument_details.tsv')
)
write_tsv(instrument_details_df, instrument_details_filename)
print(paste("Instrument details saved to:", instrument_details_filename))

#MendelianRandomization Addendum - computing an approximation of the first-stage F-statistic from genetic variants -> exposure
mr_object <- MendelianRandomization::mr_input(
  bx = harmonised_data$beta.exposure,
  bxse = harmonised_data$se.exposure,
  by = harmonised_data$beta.outcome,
  byse = harmonised_data$se.outcome,
  exposure = exposure_name,
  outcome = outcome_name,
  snps = harmonised_data$SNP
)
MRpackage_ivw_results <- MendelianRandomization::mr_ivw(
  mr_object,
  model = 'default'
)

#Addendum for MRPRESSO - horizontal pleiotropy evaluation
mrpresso_results <- tryCatch(
  {
    # This is the expression to 'try'
    mr_presso(
      data = harmonised_data,
      BetaOutcome = 'beta.outcome',
      BetaExposure = 'beta.exposure',
      SdOutcome = 'se.outcome',
      SdExposure = 'se.exposure',
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE
    )
  },
  error = function(e) {
    # This function is called ONLY if an error occurs
    message(
      'MR-PRESSO failed - likely due to insufficient instrumental variables.'
    )

    # Optional: print the original error message from R for debugging
    message('Original R error message:')
    message(e)

    # Return a specific value (like NULL) on failure
    return(NULL)
  }
)

#Output plots and data --------------------------------------------------------------------------------
#Write out a table of results
if (
  isFALSE(dir.exists(str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/'
  )))
) {
  dir.create(str_c(output_path, exposure_name, '_to_', outcome_name, '/')) #Create the directory if it doesn't exist already
}

#This outputs the print output as a txt file for the MendelianRandomization
sink(
  file = str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_MRpackage_IVW_results.txt'
  )
)
print(MRpackage_ivw_results)
sink(file = NULL)

#This outputs the print output as a txt file for the MRPresso Pleiotropy Test
if (!is.null(mrpresso_results)) {
  mrpresso_fileout <- str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_MRPRESSO_test.rds'
  )

  mrpresso_results$`MR-PRESSO results`[[
    'Outlier Test'
  ]] <- mrpresso_results$`MR-PRESSO results`[['Outlier Test']] %>%
    bind_cols('SNP' = harmonised_data$SNP)

  saveRDS(mrpresso_results, mrpresso_fileout)
}

#Save
ggsave(
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_scatter_mainline.png'
  ),
  scatter[[1]],
  dpi = 600
)
ggsave(
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_scatter_sensitivity.png'
  ),
  scatter2[[1]],
  dpi = 600
)
ggsave(
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_forest.png'
  ),
  forest[[1]],
  dpi = 600
)
ggsave(
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_looforest.png'
  ),
  loo_forest[[1]],
  dpi = 600
)

write.table(
  mainline_mr_results,
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_results_mainline.tsv'
  ),
  sep = '\t'
)
write.table(
  sensitivity_mr_results,
  str_c(
    output_path,
    exposure_name,
    '_to_',
    outcome_name,
    '/',
    exposure_name,
    '_mr_results_sensitivity.tsv'
  ),
  sep = '\t'
)
if (!is.null(or_result$result)) {
  write.table(
    or_result$result,
    str_c(
      output_path,
      exposure_name,
      '_to_',
      outcome_name,
      '/',
      exposure_name,
      '_mr_results_sensitivity_oddsratio.tsv'
    ),
    sep = '\t'
  )
} else {
  print('Odds ratio function failed')
}

#Print out a HTML report
print('Constructing HTML Report')
safe_mrreport <- safely(mr_report)
safe_mrreport(
  harmonised_data,
  output_path = str_c(output_path, exposure_name, '_to_', outcome_name, '/')
) #Generate html report
