#Script to run pleiotropy checker at scale via Rscript using Open Target Genetics'
## This is meant to be run on a internet-access node!
#Author: Jonathan Chan
#Date: 2025-05-08

library(tidyverse)
library(ggrepel)
library(otargen)
library(qvalue)
library(LDlinkR)
library(writexl)

theme_set(theme_classic())

#IMPORTER----------------------------------------------------------------------------------
#GWAS Catalog Format or Excel format for CVD HugeAmp SNPs
#Ultimate output = list of SNPs in rsID format

importer <- function(path) {
  if (str_detect(path, '\\_associations.tsv')) {
    #GWAS Catalog Format
    input_assocs <- read_tsv(path) %>%
      select(
        SNPS,
        `STRONGEST SNP-RISK ALLELE`,
        `RISK ALLELE FREQUENCY`,
        `P-VALUE`,
        `OR or BETA`
      )
    colnames(input_assocs) <- c('rsID', 'EA', 'EAF', 'pval', 'OR_beta')

    input_assocs <- input_assocs %>%
      mutate(EA = str_match(EA, '-([A-Za-z]+)')[, 2])
  } else if (str_detect(path, '_instrument_details\\.tsv')) {
    input_assocs <- read_tsv(path) %>%
      select(
        SNP,
        effect_allele.exposure,
        eaf.exposure,
        beta.exposure,
        pval.exposure
      )
    colnames(input_assocs) <- c('rsID', 'EA', 'EAF', 'OR_beta', 'pval')
  } else if (str_detect(path, '\\.xlsx')) {
    #CVD HugeAmp format
    input_assocs <- readxl::read_excel(path)
    colnames(input_assocs) <- c('rsID')
  }

  return(input_assocs)
}

#PLEIOTROPY_CHECKER-----------------------------------------------------------------
otg_caller_processor <- function(input_snp) {
  # Call OTG
  safe_pheWAS <- safely(pheWAS)
  otg_phewas_data <- safe_pheWAS(input_snp)

  # Extract out important information e.g p-value; beta
  if (is.null(otg_phewas_data$result)) {
    # i.e. no information returned
    print(str_c('No OTG PheWAS data for SNP ', input_snp))
    return(NULL) # Return NULL or an empty tibble if no data
  } else {
    processed_phewas_data <- as_tibble(otg_phewas_data$result) %>%
      mutate(SNP = input_snp)

    # Define all expected columns and their corresponding names in the input data
    # The order here defines the output column order
    expected_columns <- c(
      "SNP" = NA, # This is added by mutate, so just a placeholder for structure
      "StudyID" = "study.studyId",
      "Trait" = "study.traitReported",
      "pvalue" = "pval",
      "Beta" = "beta",
      "PMID" = "study.pmid",
      "StudySource" = "study.source",
      "StudyDate" = "study.pubDate",
      "N.Overall" = "nTotal"
    )

    # Create an empty tibble with all expected columns
    # This ensures all columns are present, even if not in the original data
    final_tibble <- tibble(.rows = nrow(processed_phewas_data))
    for (col_name in names(expected_columns)) {
      if (col_name == "SNP") {
        final_tibble[[col_name]] <- processed_phewas_data$SNP
      } else if (expected_columns[col_name] %in% names(processed_phewas_data)) {
        final_tibble[[col_name]] <- processed_phewas_data[[expected_columns[
          col_name
        ]]]
      } else {
        final_tibble[[col_name]] <- NA_real_ # Use NA_real_ for numeric, NA_character_ for character, etc.
        # Or simply NA if you don't care about type upfront
      }
    }

    # Ensure the correct column types if needed (optional, but good practice)
    # You might want to refine this based on the expected type of each column
    final_tibble <- final_tibble %>%
      mutate(
        SNP = as.character(SNP),
        StudyID = as.character(StudyID),
        Trait = as.character(Trait),
        pvalue = as.numeric(pvalue),
        Beta = as.numeric(Beta),
        PMID = as.character(PMID), # PMIDs can be numeric or character based on source
        StudySource = as.character(StudySource),
        StudyDate = as.character(StudyDate), # Or as.Date if format is consistent
        N.Overall = as.numeric(N.Overall)
      )

    return(final_tibble)
  }
}

#This applies multiple testing correction to filter for only MTC-significant pleiotropic associations.

fdr_corrector <- function(processed_data_tb, fdr_cutoff = 0.05) {
  #Generate a FDR for each study's association taking into account that OTG only reports studies with a pvalue < 0.005 so need to use qvalue_truncp() as in qvalue R package as in Harper et al, 2021 Supp. Table 25

  #This outputs the adjusted p-value for each study within a SNP. The input_tb is on a per-SNP basis.
  safe_qvalue_truncp <- safely(qvalue_truncp)

  qvals <- safe_qvalue_truncp(
    processed_data_tb$pvalue,
    fdr.level = fdr_cutoff,
    pfdr = T
  )
  qvals <- qvals$result

  if (is_null(qvals)) {
    output <- processed_data_tb %>%
      mutate(qvalue = NA)
  } else {
    output <- bind_cols(processed_data_tb, qvalue = qvals$qvalues)
  }

  return(output)
}

scatter_plotter <- function(
  processed_data_tb,
  pheno,
  output_path,
  qval_cutoff = 0.05,
  priority_cutoff = 50
) {
  theme_set(theme_classic())

  snp <- unique(processed_data_tb$SNP)
  processed_data_tb2 <- processed_data_tb %>% filter(qvalue < qval_cutoff) #This filters the PheWAS results for only those where the SNP is associated at GWS with the trait.

  if (length(unique(processed_data_tb2$Trait)) >= priority_cutoff) {
    #This defines a priority file prefix if there is > n significant signals PRIOR to MTC across phenotypes btw
    prefix <- 'priority_'
  } else {
    prefix <- ''
  }

  if (nrow(processed_data_tb2) > 0) {
    out <- ggplot(processed_data_tb2, aes(y = Trait, x = -log10(qvalue))) +
      geom_point(aes(col = StudySource)) +
      ylab('Trait') +
      xlab('-log10(FDR-adjusted p-value)') +
      labs(title = str_c('PheWAS from OTG for SNP ', snp))

    print(out)
    ggsave(
      str_c(
        output_path,
        '/2_phewas_plots/',
        pheno,
        '/',
        prefix,
        snp,
        '_phewas_scatter.png'
      ),
      out,
      dpi = 600,
      create.dir = T,
      width = 15,
      height = 12
    )
  }
}

allsnp_merger_writer <- function(
  fdr_corrected_data,
  pheno,
  output_path,
  fdr_cutoff = 0.05
) {
  output <- fdr_corrected_data %>% bind_rows() %>% filter(qvalue < fdr_cutoff) #This requires the adjusted p-value to be at least 0.01

  write_tsv(output, str_c(output_path, pheno, '_phewas_FDRfiltered.tsv'))

  return(output)
}

pairwise_LD_checker <- function(
  input_SNPs,
  pheno,
  output_path,
  token,
  population = 'GBR'
) {
  #Note this API token means don't push this to Git

  all_possible_combinations <- expand_grid(input_SNPs, input_SNPs)
  colnames(all_possible_combinations) <- c('var1', 'var2')
  all_possible_combinations <- filter(all_possible_combinations, var1 != var2)

  safe_LDpair <- safely(LDpair)
  pairwise_ld_results <- map2(
    all_possible_combinations$var1,
    all_possible_combinations$var2,
    ~ safe_LDpair(
      .x,
      .y,
      pop = population,
      token,
      genome_build = 'grch38_high_coverage'
    )[['result']]
  ) %>%
    bind_rows() %>%
    filter(!is.na(p_val) & p_val < 0.05)

  write_tsv(pairwise_ld_results, str_c(output_path, pheno, '_pairwiseLD.tsv'))

  return(pairwise_ld_results)
}

#This is the main pleiotropy function
main2 <- function(path, pheno, output_path, pairwise_LD_check = F) {
  print(str_c('Importing in ', path))
  input_assocs <- importer(path)

  if (
    file.exists(str_c(output_path, '/1_OTG_RDS/', pheno, '_phewas_results.rds'))
  ) {
    processed_data <- read_rds(str_c(
      output_path,
      '/1_OTG_RDS/',
      pheno,
      '_phewas_results.rds'
    ))
  } else {
    processed_data <- map(input_assocs$rsID, ~ otg_caller_processor(.))
    processed_data <- processed_data[map_lgl(processed_data, ~ !is.null(.))]
    write_rds(
      processed_data,
      str_c(output_path, '/1_OTG_RDS/', pheno, '_phewas_results.rds')
    )
  }

  fdr_corrected_data <- map(processed_data, ~ fdr_corrector(.)) #This performs MTC for the pleiotropy analysis
  walk2(fdr_corrected_data, pheno, ~ scatter_plotter(.x, .y, output_path))

  fdrsig_signals <- allsnp_merger_writer(fdr_corrected_data, pheno, output_path)

  #Output a summary table of the number of unique FDR-significant traits per SNP
  summary_table <- fdrsig_signals %>%
    group_by(SNP) %>%
    summarize(NumberOfUniqueTraits = length(unique(Trait))) %>%
    ungroup() %>%
    arrange(desc(NumberOfUniqueTraits))

  write_tsv(
    summary_table,
    str_c(output_path, pheno, '_phewas_FDRfiltered_summaryN.tsv')
  )

  if (isTRUE(pairwise_LD_check)) {
    token <- readRDS('ldlinkR_API_token.rds')
    pairwise_LD_checker(input_assocs$rsID, pheno, output_path, token)
  }
}

snp_pleio_filter <- function(
  og_path,
  pheno,
  output_path,
  red_flags = '',
  pheno2 = '',
  n_cutoff = 50
) {
  input_tb <- read_tsv(str_c(output_path, pheno, '_phewas_FDRfiltered.tsv'))

  #These are the SNPs which have pleiotropy with a red flag trait which is known to be causally associated with the outcome
  if (!identical(red_flags, '')) {
    redflag_snps <- input_tb %>%
      group_by(SNP) %>%
      filter(str_detect(Trait, str_c(red_flags, collapse = '|'))) %>%
      select(SNP) %>%
      unique()
  } else {
    redflag_snps <- input_tb %>%
      filter(is.infinite(1)) #i.e empty tibble
  }

  #These are the SNPs which have excess pleiotropy with many other traits
  n_failed_snps <- input_tb %>%
    group_by(SNP) %>%
    select(SNP, Trait) %>%
    mutate(NumberOfUniqueTraits = length(unique(Trait))) %>%
    ungroup() %>%
    filter(NumberOfUniqueTraits > n_cutoff) %>%
    select(SNP) %>%
    unique()

  #Output the SNPs which fail the pleiotropy checks
  output_snps <- unique(c(redflag_snps$SNP, n_failed_snps$SNP))

  og_snps <- importer(og_path)

  print(str_c(
    'For phenotype ',
    pheno,
    ' of ',
    length(unique(og_snps$rsID)),
    ' instruments, ',
    'total of ',
    length(output_snps),
    ' instruments excluded with ',
    nrow(n_failed_snps),
    ' instruments excluded due to excessive pleiotropy with >',
    n_cutoff,
    ' traits and ',
    nrow(redflag_snps),
    ' instruments excluded due to association with traits linked to the oucome e.g. ',
    str_c(red_flags, collapse = '|')
  ))

  #Write out the tsv of rsIDs of SNPs which pass pleiotropy checks as well
  pleio_pass <- og_snps %>%
    filter(!rsID %in% output_snps) %>%
    select(rsid = rsID) %>%
    unique()

  if (pheno != 'HCM') {
    write_xlsx(
      pleio_pass,
      str_c(
        '../1_data/gwas_associations/1_pleio_pruned/',
        pheno,
        '_gwas_associations.xlsx'
      )
    )
  } else {
    write_xlsx(
      pleio_pass,
      str_c(
        '../1_data/gwas_associations/1_pleio_pruned/',
        pheno,
        '_to_',
        pheno2,
        '_gwas_associations.xlsx'
      )
    )
  }

  #Also write out TSV of failed SNPs
  write_tsv(
    tibble(SNP = output_snps) %>% unique(),
    str_c(output_path, pheno, '_pleiofailedsnps.tsv')
  )
}
# RUN -----------------------------------------------------------------------------------
instrument_files <- c(
  '../1_data/gwas_associations/tadros25_hcm_nonMTAG_associations.xlsx',
  '../1_data/gwas_associations/roselli25_af_associations.xlsx',
  '../1_data/gwas_associations/keaton24_dbp_associations.tsv',
  '../1_data/gwas_associations/sidorenko24_bmi_associations.tsv',
  '../1_data/gwas_associations/sakaue21_t2d_associations.tsv',
  '../1_data/gwas_associations/liu18_SmokingInitiation_associations.tsv',
  '../1_data/gwas_associations/nauffal24_t1time_associations.xlsx',
  '../1_data/gwas_associations/thanaj22_maxLAVi_associations.tsv',
  '../1_data/gwas_associations/thanaj22_long_dsr_associations.tsv',
  '../1_data/gwas_associations/thanaj22_radial_dsr_associations.tsv'
)
phenotypes <- c(
  'HCM',
  'AF',
  'DBP',
  'BMI',
  'T2D',
  'SmokingInitiation',
  'T1time',
  'maxLAVi',
  'longDSR',
  'radialDSR'
)

walk2(
  instrument_files,
  phenotypes,
  ~ main2(.x, .y, '../3_output/1_MR_QC/pleiotropy_check/')
)

hcmr_instrument_files <- list.files(
  '../1_data/gwas_associations/',
  pattern = 'hcmr[^\\.]+\\.xlsx',
  full.names = T
)
hcmr_phenotypes <- str_match(
  basename(hcmr_instrument_files),
  '([^\\.]+)_associations\\.xlsx'
)[, 2]

ukb_instrument_files <- list.files(
  '../1_data/gwas_associations/',
  pattern = 'ukb[^\\.]+\\.xlsx',
  full.names = T
)
ukb_phenotypes <- str_match(
  basename(ukb_instrument_files),
  '([^\\.]+)_associations\\.xlsx'
)[, 2]

walk2(
  hcmr_instrument_files,
  hcmr_phenotypes,
  ~ main2(.x, .y, '../3_output/1_MR_QC/pleiotropy_check/')
)
walk2(
  ukb_instrument_files,
  ukb_phenotypes,
  ~ main2(.x, .y, '../3_output/1_MR_QC/pleiotropy_check/')
)

# PLEIOTROPY FILTER ---------------------------------------------------------------

#DEFINE THE RED FLAG TRAITS
#For the X -> HCM
hcm_outcome_redflags <- c(
  'Body mass index',
  'Pulse pressure',
  'Systolic blood pressure',
  'Diastolic blood pressure',
  'Hypertension',
  'Hypertensive diseases'
)
#There are obviously exceptions, when the exposure is BMI or DBP
bmi_exception <- c(
  'Pulse pressure',
  'Systolic blood pressure',
  'Diastolic blood pressure',
  'Hypertension',
  'Hypertensive diseases'
)
dbp_exception <- c('Body mass index')

## OVERALL RED FLAGS
redflags <- c(
  list('', hcm_outcome_redflags, dbp_exception, bmi_exception),
  rep(list(hcm_outcome_redflags), length(phenotypes) - 4)
) #i.e no red flag traits for HCM instruments
hcmr_redflags <- rep(list(hcm_outcome_redflags), length(hcmr_phenotypes))
ukb_redflags <- rep(list(hcm_outcome_redflags), length(ukb_phenotypes))
rm(bmi_exception, dbp_exception)

#Run for public GWAS summstats
pwalk(
  list(instrument_files, phenotypes, redflags),
  ~ snp_pleio_filter(..1, ..2, '../3_output/1_MR_QC/pleiotropy_check/', ..3)
)

#Run for internal HCMR and UKB GWAS summstats
pwalk(
  list(hcmr_instrument_files, hcmr_phenotypes, hcmr_redflags),
  ~ snp_pleio_filter(..1, ..2, '../3_output/1_MR_QC/pleiotropy_check/', ..3)
)
pwalk(
  list(ukb_instrument_files, ukb_phenotypes, ukb_redflags),
  ~ snp_pleio_filter(..1, ..2, '../3_output/1_MR_QC/pleiotropy_check/', ..3)
)


# MR YAML WRITER----------------------------------------------------------------
#' Create a YAML configuration file for Mendelian Randomization analysis.

create_mr_yaml_config <- function(
  non_hcm_phenotype,
  non_hcm_exposure_gwas_path,
  non_hcm_instruments_path,
  output_yaml_file_path,
  hcm_instruments_path = '../1_data/gwas_associations/1_pleio_pruned/HCM_to__gwas_associations.xlsx'
) {
  # Required package for YAML writing
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "Package 'yaml' is required but not installed. Please install it using install.packages('yaml')."
    )
  }
  # No need to explicitly call library(yaml) if using ::, but good for clarity if used frequently
  # library(yaml)

  # --- Define static values based on the provided YAML structure ---
  # These values are fixed according to your example
  static_output_path_val <- '../3_output/2_MR/'
  static_ld_eur_bed_file_val <- '../../hcmr_ukbb/popgen/5_MR/pipeline_files/1kg_v3_ld/EUR'
  static_plink_binary_path_val <- '/well/PROCARDIS/jchan/bin/plink'

  # Details for d1 (Direction 1: non-HCM phenotype -> HCM)
  static_d1_outcome_phenotype_val <- 'HCM'
  static_d1_outcome_path_val <- 'tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'
  static_d1_ld_clump_val <- 'FALSE' # R logical FALSE becomes 'false' in YAML output

  # Details for d2 (Direction 2: HCM -> non-HCM phenotype)
  static_d2_exposure_phenotype_val <- 'HCM'
  static_d2_exposure_path_val <- 'tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'

  # Construct the path for d2_instruments dynamically
  # This path is for instruments when HCM is the exposure and non-HCM phenotype is the outcome
  # Example: '../1_data/gwas_associations/1_pleio_pruned/HCM_to_BMI_gwas_associations.xlsx'
  d2_instruments_path_val <- hcm_instruments_path
  static_d2_ld_clump_val <- 'FALSE' # R logical FALSE becomes 'false' in YAML output

  # --- Construct the list that represents the YAML structure ---
  config_data_list <- list(
    output_path = static_output_path_val,
    ld_eur_bed_file = static_ld_eur_bed_file_val,
    plink_binary_path = static_plink_binary_path_val,

    # Direction 1: non-HCM phenotype (exposure) -> HCM (outcome)
    d1_exposure = non_hcm_phenotype,
    d1_outcome = static_d1_outcome_phenotype_val,
    d1_exposure_path = non_hcm_exposure_gwas_path, # User-provided path for non-HCM GWAS
    d1_outcome_path = static_d1_outcome_path_val, # Fixed path for HCM GWAS
    d1_instruments = non_hcm_instruments_path, # User-provided path for non-HCM instruments
    d1_ld_clump = static_d1_ld_clump_val,

    # Direction 2: HCM (exposure) -> non-HCM phenotype (outcome)
    d2_exposure = static_d2_exposure_phenotype_val,
    d2_outcome = non_hcm_phenotype, # Non-HCM phenotype is the outcome here
    d2_exposure_path = static_d2_exposure_path_val, # Fixed path for HCM GWAS (as exposure)
    d2_outcome_path = non_hcm_exposure_gwas_path, # User-provided path for non-HCM GWAS (as outcome)
    d2_instruments = d2_instruments_path_val, # Dynamically constructed path
    d2_ld_clump = static_d2_ld_clump_val
  )

  # --- Write the list to a YAML file ---
  tryCatch(
    {
      yaml::write_yaml(config_data_list, file = output_yaml_file_path)
      message(paste0(
        "Successfully created YAML configuration file at: ",
        output_yaml_file_path
      ))
      invisible(TRUE) # Return TRUE invisibly on success
    },
    error = function(e) {
      # Provide a more informative error message
      stop(paste0(
        "Failed to write YAML file to '",
        output_yaml_file_path,
        "'. Error: ",
        e$message
      ))
    }
  )
}

hcmr_instrument_pleiopruned <- str_c(
  dirname(hcmr_instrument_files),
  '/1_pleio_pruned/',
  hcmr_phenotypes,
  '_gwas_associations.xlsx'
)
ukb_instrument_pleiopruned <- str_c(
  dirname(ukb_instrument_files),
  '/1_pleio_pruned/',
  ukb_phenotypes,
  '_gwas_associations.xlsx'
)

hcmr_gwas_summstats <- str_c(
  hcmr_phenotypes,
  '_gwas_ssf/final/',
  hcmr_phenotypes,
  '_gwas_ssf.h.tsv.gz'
)
ukb_gwas_summstats <- str_c(
  ukb_phenotypes,
  '_gwas_ssf/final/',
  ukb_phenotypes,
  '_gwas_ssf.h.tsv.gz'
)

hcmr_yaml_outputpath <- str_c(
  '../1_data/input_configs/',
  hcmr_phenotypes,
  '_and_HCM.yaml'
)
ukb_yaml_outputpath <- str_c(
  '../1_data/input_configs/',
  ukb_phenotypes,
  '_and_HCM.yaml'
)


pwalk(
  list(
    hcmr_phenotypes,
    hcmr_gwas_summstats,
    hcmr_instrument_pleiopruned,
    hcmr_yaml_outputpath
  ),
  ~ create_mr_yaml_config(..1, ..2, ..3, ..4)
)

pwalk(
  list(
    ukb_phenotypes,
    ukb_gwas_summstats,
    ukb_instrument_pleiopruned,
    ukb_yaml_outputpath
  ),
  ~ create_mr_yaml_config(..1, ..2, ..3, ..4)
)
