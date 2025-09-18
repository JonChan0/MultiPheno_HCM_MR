'''
This outlines the Snakemake pipeline to harmonise all input GWAS summary statistics using EBI GWAS Catalog's GWAS Summstats Harmoniser.
This version includes a conditional rule to reformat specific GWAS outputs instead of harmonising them.
Author: Jonathan Chan
Date: 2025-02-10
'''

import os
import pandas as pd
import hashlib
import yaml

configfile: 'config/config_gwas_summstats_harmoniser.yaml'

# --- Input File Partitioning ---
# This section separates input files into two groups based on their path.

HCMR_PATH = '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/hcmr/'
UKB_PATH = '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/ukb/'
all_input_paths = config['input_tsv_gz_files']

# HCMR GWAS files to be reformatted which have the format of phenotype_manhattan_rsid2.tsv
HCMR_paths = [p for p in all_input_paths if HCMR_PATH in p]
HCMR_folders = [os.path.dirname(p) for p in HCMR_paths]
#Grab the basename without the full directory path, only the filename without the extension but also without _manhattan_rsid2.tsv
HCMR_basenames = [os.path.splitext(os.path.splitext(os.path.basename(p))[0])[0].replace('_manhattan_rsid2', '') for p in HCMR_paths]

UKB_paths = [p for p in all_input_paths if UKB_PATH in p]
UKB_folders = [os.path.dirname(p) for p in UKB_paths]
#Grab the basename without the full directory path, only the filename without the extension but also without _manhattan_rsid2.tsv
UKB_basenames = [os.path.splitext(os.path.splitext(os.path.basename(p))[0])[0].replace('_manhattan_rsid', '') for p in UKB_paths]

# Standard files to be harmonised which have the format of phenotype.tsv.gz
standard_paths = [p for p in all_input_paths if HCMR_PATH not in p and UKB_PATH not in p]
standard_folders = [os.path.dirname(p) for p in standard_paths]
standard_basenames = [os.path.splitext(os.path.splitext(os.path.basename(p))[0])[0] for p in standard_paths]

# Create a complete list of all basenames and a map to find their original folder
all_basenames = standard_basenames + HCMR_basenames + UKB_basenames
all_folders = standard_folders + HCMR_folders + UKB_folders
basename_to_folder_map = dict(zip(all_basenames, all_folders))

hcmr_only_output = expand('{folder}/{basename}_gwas_ssf.tsv.gz', zip, folder=HCMR_folders, basename=HCMR_basenames)
ukb_only_output = expand('{folder}/{basename}_gwas_ssf.tsv.gz', zip, folder=UKB_folders, basename=UKB_basenames)
common_output = [
    expand(config['base_output_folder']+'{basename}_gwas_ssf/final/{basename}_gwas_ssf.h.tsv.gz', basename=all_basenames),
    config['desired_output_folder']+'harmonisation_summary.tsv',
    config['desired_output_folder']+'max_n_summary.tsv'
]

rule all:
    input: 
        hcmr_only_output + ukb_only_output + common_output

# --- Specialised Rule for OG Data ---
# Rule to reformat 'manhattan_rsid2.tsv' from HCMR or 'manhattan_rsid.tsv' from UKB into the GWAS-SSF format.
# This rule runs ONLY for inputs which have the HCMR_PATH or UKB_PATH in their path.
rule reformat_hcmr_gwas:
    input:
        raw_gwas = lambda wildcards: (
            f"{wildcards.folder}/{wildcards.basename}_manhattan_rsid2.tsv"
            if "hcmr" in wildcards.basename.lower()
            else f"{wildcards.folder}/{wildcards.basename}_manhattan_rsid.tsv"
        )
    output:
        reformatted_gwas_gz = "{folder}/{basename}_gwas_ssf.tsv.gz",
        meta_yaml = "{folder}/{basename}_gwas_ssf.tsv.gz-meta.yaml"
    resources:
        mem_mb = 16000
    shell:
        """
        # Use a python script for robust reformatting and YAML creation
        python -c "
import pandas as pd
import hashlib
import os

# To escape braces for Snakemake, you must double them: {{ }}
col_map = {{
    'rsid': 'variant_id', 'chromosome': 'chromosome', 'position': 'base_pair_location',
    'allele_B': 'effect_allele', 'allele_A': 'other_allele', 'eaf': 'effect_allele_frequency',
    'beta': 'beta', 'beta_se': 'standard_error', 'pval': 'p_value', 'all_total': 'n_total'
}}

# Read and reformat data
df = pd.read_csv('{input.raw_gwas}', sep='\\s+') # Using whitespace separator
df_reformatted = df.rename(columns=col_map)[list(col_map.values())]
df_reformatted.to_csv('{output.reformatted_gwas_gz}', sep='\\t', index=False, compression='gzip', na_rep='NA')

# Calculate MD5 checksum
with open('{output.reformatted_gwas_gz}', 'rb') as f:
    md5_hash = hashlib.md5(f.read()).hexdigest()

# Create YAML metadata file
# Here, we escape the f-string braces with {{ }} so Snakemake ignores them.
# The inner {output.reformatted_gwas_gz} is left as-is so Snakemake can format it.
yaml_content = f'''# Study meta-data
date_metadata_last_modified: 2023-02-09

# Genotyping Information
genome_assembly: GRCh37
coordinate_system: 1-based
genotyping_technology:
  - Genome-wide genotyping array

# Summary Statistic information
data_file_name: {{os.path.basename(r'{output.reformatted_gwas_gz}')}}
file_type: GWAS-SSF v0.1
data_file_md5sum: {{md5_hash}}

# Harmonization status
is_harmonised: false
is_sorted: false
'''
with open('{output.meta_yaml}', 'w') as f:
    f.write(yaml_content)
"
        """

#Rule 1 = Standardise the headings and columns via gwas-ssf by GWAS Summstats Tools
# rule gwas_ssf:
#     input: "{folder}/{basename}.tsv.gz"
#     output:
#         gwas_sumstats_json = temp("{folder}/{basename}.json"),
#         standardised_gwas_summstats_tsv = temp("{folder}/{basename}_gwas_ssf.tsv"),
#         standardised_gwas_summstats_tsv_gz = "{folder}/{basename}_gwas_ssf.tsv.gz"
#     conda: 'gms'
#     resources:
#         mem_mb=8000
#     shell:
#         '''
#         gwas-ssf format {input} --generate_config --config_out {output.gwas_sumstats_json}
#         gwas-ssf format {input} --apply_config --config_in {output.gwas_sumstats_json} -o {output.standardised_gwas_summstats_tsv}
#         gzip {output.standardised_gwas_summstats_tsv} -c > {output.standardised_gwas_summstats_tsv_gz}
#         '''

rule gwas_summstats_harmoniser:
    input:
        # The lambda function uses the map to find the correct input path for any given basename.
        gwas_summstats=lambda wildcards: f"{basename_to_folder_map[wildcards.basename]}/{wildcards.basename}_gwas_ssf.tsv.gz"
    output:
        harmonised_output_folder = directory(config['base_output_folder']+'{basename}_gwas_ssf/'),
        harmonised_gwas_summstats = config['base_output_folder']+'{basename}_gwas_ssf/final/{basename}_gwas_ssf.h.tsv.gz',
        output_logfiles = config['base_output_folder']+'{basename}_gwas_ssf/final/{basename}_gwas_ssf.running.log'
    resources:
        mem_mb = 40000
    params:
        ref_folder = config['gwas_summstats_harmoniser_ref_folder'],
        output_folder = config['desired_output_folder']
    shell:
        '''
        module load Nextflow/24.04.2

        nextflow run EBISPOT/gwas-sumstats-harmoniser -r v1.1.10 \
        --ref {params.ref_folder} \
        --file {input.gwas_summstats} \
        --harm \
        --chromlist 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
        --profile executor,conda \
        --to_build '38'

        echo Successfully completed harmonisation of {input.gwas_summstats}
        '''

# Rule 3 = Parse all '.running.log' files to generate a single summary report
rule summarise_harmonisation:
    input:
        # Collect log files from ALL harmonisation jobs.
        lHCMR_files = expand(config['base_output_folder']+'{basename}_gwas_ssf/final/{basename}_gwas_ssf.running.log', basename=all_basenames)
    output:
        summary_file = config['desired_output_folder']+'harmonisation_summary.tsv'
    run:
        summary_data = []
        for lHCMR_file in input.lHCMR_files:
            basename = os.path.basename(lHCMR_file).split('_gwas_ssf')[0]
            with open(lHCMR_file, 'r') as f:
                lines = f.readlines()
                success = any("Result\tSUCCESS_HARMONIZATION" in line for line in lines)
                if success:
                    percentage_line = [line for line in lines if 'sites successfully harmonised' in line]
                    clean_percentage = percentage_line[0].strip() if percentage_line else "N/A"
                    summary_data.append(f"{basename}\tYes\t{clean_percentage}")
                else:
                    summary_data.append(f"{basename}\tNo\t(N/A)")
        
        with open(output.summary_file, 'w') as f:
            f.write("Phenotype\tHarmonised\tPercentage Harmonised\n")
            for line in sorted(summary_data): # Sort for consistent output
                f.write(line + "\n")

# Rule 4 = Extract maximum sample size (N) using a scatter-gather bash approach.

# Rule 4a = (Scatter) Extract max N for each file individually using awk.
# This rule processes one GWAS file and outputs one single-row temporary file.
rule extract_max_n_per_file:
    input:
        harmonised_file=config['base_output_folder']+'{basename}_gwas_ssf/final/{basename}_gwas_ssf.h.tsv.gz'
    output:
        # Create a temporary file for each input to hold its max N value
        max_n_row=config['desired_output_folder']+'temp_max_n/{basename}.max_n.tsv'
    shell:
        r"""
        # This command unzips the file and pipes it to awk.
        # 1. BEGIN: Set the field separator to a tab, initialize max=0 and column index=-1.
        # 2. NR==1: On the header row, loop through columns to find 'n' or 'n_total' (case-insensitive) and store its index.
        # 3. NR>1: For all data rows, if the sample size column contains a number larger than the current max, update max.
        # 4. END: After processing all lines, print the trait name (passed from wildcards) and the final max value.
        # The output is a single row, e.g., "trait_name\t123456".
        zcat {input.harmonised_file} | awk -v trait="{wildcards.basename}" \
        'BEGIN {{ FS="\t"; max=0; col_idx=-1 }} \
        NR==1 {{ for(i=1; i<=NF; i++) {{ if (tolower($i) == "n" || tolower($i) == "n_total") {{ col_idx=i; break }} }} }} \
        NR>1 {{ if (col_idx != -1 && $col_idx ~ /^[0-9.]+$/ && $col_idx > max) {{ max=$col_idx }} }} \
        END {{ if (col_idx == -1) {{ print trait"\tNA" }} else {{ printf "%s\t%d\n", trait, max }} }}' \
        > {output.max_n_row}
        """

# Rule 4b = (Gather) Aggregate all individual max N files into the final summary file.
# This rule collects all the temporary files and combines them.
rule aggregate_max_n:
    input:
        # Collect all the temporary files generated by the rule above
        max_n_files=expand(config['desired_output_folder']+'temp_max_n/{basename}.max_n.tsv', basename=all_basenames)
    output:
        max_n_summary=config['desired_output_folder']+'max_n_summary.tsv'
    shell:
        """
        # First, write the header to the final output file
        echo -e "trait\tmax_N" > {output.max_n_summary}
        # Then, concatenate all the single-row temporary files and append them to the final output
        cat {input.max_n_files} >> {output.max_n_summary}
        """