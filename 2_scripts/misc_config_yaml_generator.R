#Script to generate the config.yamls for 2_bidirectionalMR.smk 
#Author: Jonathan Chan
#Date: 2025-02-20

library(tidyverse)


#Define the fixed variables here
output_path <- '../3_output/2_MR/'
ld_file <- '../../hcmr_ukbb/popgen/5_MR/pipeline_files/1kg_v3_ld/EUR'
plink_file <- '/well/PROCARDIS/jchan/bin/plink'

assoc_basepath <- '../1_data/gwas_associations/'

#Grab out the phenotype names from the folders in the working directory
phenos <- list.dirs('./', recursive=F)
phenos <- phenos[str_detect(phenos, 'gwas_ssf')]
phenos <- phenos[!str_detect(phenos,'hcm')] #Drop HCM
phenos2 <- str_match(phenos, '_([A-Za-z_\\d]+)_gwas_ssf')[,2]


#Define the fixed HCM related 
d1_outcome <- 'HCM'
d2_exposure <- 'HCM'
d1_outcome_path <- d2_exposure_path <- 'tadros24_hcm_gwas_ssf/final/tadros24_hcm_gwas_ssf.h.tsv.gz'
d2_instruments <- "rs1048302,rs11687178,rs3845778,rs2540277,rs6747402,rs7612736,rs4894803,rs2191446,rs11748963,rs3176326,rs12210733,rs66520020,rs7824244,rs35006907,rs2645210,rs2177843,rs11196085,rs17617337,rs12270374,rs182427065,rs7487962,rs41306688,rs113907726,rs8006225,rs8033459,rs28768976,rs7210446,rs2644262,rs6566955,rs12460541,rs62222424,rs5760054"
d2_ld_clump <- F

print(str_c('Number of HCM instruments = ', length(str_split(d2_instruments,',')[[1]])))

#Define the variable phenotype related
d1_exposure_path <- str_c(str_match(phenos, '([A-Za-z_\\d]+)_gwas_ssf')[,2],'_gwas_ssf/final/',str_match(phenos, '([A-Za-z_\\d]+)_gwas_ssf')[,2],'_gwas_ssf.h.tsv.gz')
d1_exposure <-phenos2
d2_outcome <- phenos2

assocs <- str_c(assoc_basepath, list.files(assoc_basepath)[-str_detect(list.files(assoc_basepath), 'archive')])

variables <- tibble(d1_exposure=phenos2, d2_outcome=phenos2, d1_exposure_path = d1_exposure_path)
variables <- variables %>%
  rowwise() %>%
  mutate(d1_instruments = ifelse(any(str_detect(assocs,str_c('_',d1_exposure,'_'))), 
                                 assocs[str_detect(assocs, str_c('_',d1_exposure,'_'))],NA)) %>%
  mutate(d1_ld_clump = ifelse(str_detect(d1_instruments, '\\.csv'), T, F)) %>% #If it detects .csv, then clump the SNPs
  mutate(d1_exposure = case_when(
    d1_exposure %in% c('crp', 'hf','af','t2d','rvef','rvedvi','rvesvi','rvsvi')~ str_to_upper(d1_exposure),
    d1_exposure == 'rafac' ~ 'RA_Frac_Change',
    d1_exposure =='t1time' ~ 'T1time',
    T~d1_exposure)) %>%
  mutate(d2_outcome = case_when(
    d2_outcome %in% c('crp', 'hf','af','t2d','rvef','rvedvi','rvesvi','rvsvi')~ str_to_upper(d1_exposure),
    d2_outcome == 'rafac' ~ 'RA_Frac_Change',
    d2_outcome == 't1time'~'T1time',
    T~d2_outcome) )

#Remove DCM row
variables <- filter(variables, d1_exposure != 'dcm')
                                 
#Add in the fixed variables
variables <- variables %>%
  mutate(output_path = output_path,
         ld_eur_bed_file = ld_file,
         plink_binary_path = plink_file, 
         d1_outcome = d1_outcome,
         d2_exposure = d2_exposure,
         d1_outcome_path = d1_outcome_path,
         d2_exposure_path = d2_exposure_path,
         d2_outcome_path = d1_exposure_path, #Because the d1_exposure == d2_outcome
         d2_instruments = d2_instruments, 
         d2_ld_clump = d2_ld_clump)

#Write out a single YAML file for each row in the tibble

yaml_output <- function(tb_row, output_yaml_path='../1_data/input_configs/'){
  
  # Create the lines for the output YAML file
  lines <- c(
    paste0("output_path: '", tb_row$output_path, "'"),
    paste0("ld_eur_bed_file: '", tb_row$ld_eur_bed_file, "'"),
    paste0("plink_binary_path: '", tb_row$plink_binary_path, "'"),
    
    paste0("d1_exposure: '", tb_row$d1_exposure, "'"),
    paste0("d1_outcome: '", tb_row$d1_outcome, "'"),
    paste0("d1_exposure_path: '", tb_row$d1_exposure_path, "'"),
    paste0("d1_outcome_path: '", tb_row$d1_outcome_path, "'"),
    paste0("d1_instruments: '", tb_row$d1_instruments, "'"),
    paste0("d1_ld_clump: '", tb_row$d1_ld_clump, "'"),
    
    paste0("d2_exposure: '", tb_row$d2_exposure, "'"),
    paste0("d2_outcome: '", tb_row$d2_outcome, "'"),
    paste0("d2_exposure_path: '", tb_row$d2_exposure_path, "'"),
    paste0("d2_outcome_path: '", tb_row$d2_outcome_path, "'"),
    paste0("d2_instruments: \"", tb_row$d2_instruments, "\""),
    paste0("d2_ld_clump: '", tb_row$d2_ld_clump, "'")
  )
  
  # Write the lines to a file
  writeLines(lines, con = str_c(output_yaml_path, tb_row$d1_exposure, '_and_', tb_row$d1_outcome, '.yaml'))
  
}

walk(split(variables, seq_len(nrow(variables))), ~yaml_output(.))



