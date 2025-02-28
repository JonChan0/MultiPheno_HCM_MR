'''
This outlines the Snakemake pipeline to harmonise all input GWAS summary statistics using EBI GWAS Catalog's GWAS Summstats Harmoniser
Author: Jonathan Chan
Date: 2025-02-10
'''

import os

configfile: 'config_gwas_summstats_harmoniser.yaml'

input_tsv_gz_paths = config['input_tsv_gz_files']
#Extract out the basename with the full directory path i.e just without the .tsv.gz extension
# input_tsv_gz_paths = [os.path.splitext(os.path.splitext(input_tsv_gz_path)[0])[0] for input_tsv_gz_path in input_tsv_gz_paths]
input_tsv_gz_folders = [os.path.dirname(x) for x in input_tsv_gz_paths]
# print(input_tsv_gz_folders)

#Extract out the basename without the full directory path, only the filename without the extension. There are two extensions to remove i.e convert ...tsv.gz to ...
input_tsv_gz_basenames = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in input_tsv_gz_paths]

rule all:
    input: 
        standardised_gwas_summstats = expand("{input_tsv_gz_folder}/{input_tsv_gz_basename}_gwas_ssf.tsv.gz", zip, input_tsv_gz_folder=input_tsv_gz_folders, input_tsv_gz_basename=input_tsv_gz_basenames),
        harmonised_gwas_summstats = expand(config['base_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/final/{input_tsv_gz_basename}_gwas_ssf.h.tsv.gz', input_tsv_gz_basename=input_tsv_gz_basenames),
        summary_file = config['desired_output_folder']+'harmonisation_summary.tsv'


#Rule 1 = Standardise the headings and columns via gwas-ssf by GWAS Summstats Tools
rule gwas_ssf:
    input: "{input_tsv_gz_folder}/{input_tsv_gz_basename}.tsv.gz"
    output:
        gwas_sumstats_json = temp("{input_tsv_gz_folder}/{input_tsv_gz_basename}.json"),
        standardised_gwas_summstats_tsv = temp("{input_tsv_gz_folder}/{input_tsv_gz_basename}_gwas_ssf.tsv"),
        standardised_gwas_summstats_tsv_gz = "{input_tsv_gz_folder}/{input_tsv_gz_basename}_gwas_ssf.tsv.gz"
    conda: 'gms'
    resources:
        mem_mb=8000
    shell:
        '''
        gwas-ssf format {input} --generate_config --config_out {output.gwas_sumstats_json}
        gwas-ssf format {input} --apply_config --config_in {output.gwas_sumstats_json} -o {output.standardised_gwas_summstats_tsv}
        gzip {output.standardised_gwas_summstats_tsv} -c > {output.standardised_gwas_summstats_tsv_gz}
        '''

#Rule 2 = Run GWAS summstats harmoniser using the above input files
rule gwas_summstats_harmoniser:
    input:
        gwas_summstats=lambda wildcards: f"{input_tsv_gz_folders[input_tsv_gz_basenames.index(wildcards.input_tsv_gz_basename)]}/{wildcards.input_tsv_gz_basename}_gwas_ssf.tsv.gz"
    output:
        harmonised_output_folder = directory(config['base_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/'),
        harmonised_gwas_summstats = config['base_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/final/{input_tsv_gz_basename}_gwas_ssf.h.tsv.gz',
        output_logfiles = config['base_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/final/{input_tsv_gz_basename}_gwas_ssf.running.log'
    # conda: 'gms'
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
        # mv {output.harmonised_output_folder} -t {params.output_folder}
        '''

#I manually move all the folders over to the desired output folder because smk not working for some reason.

#Rule 3 = Parse the '.running.log' file to check if the harmonisation has succeded correctly given the a successful harmonisation will have the string 'Result  SUCCESS_HARMONIZATION' in the file
#If this is the case, please report out the percentage of variants successfully harmonised by reporting out the string on the line prior to 'sites successfully harmonised.' e.g '97.39% ( 16263377 of 16698468 ) sites successfully harmonised.'
#Iterate over all the different phenotypes to output a single summary file with the percentage of variants successfully harmonised for each phenotype and whether or not each phenotype was successfully harmonised
#For each line in the output file, print first the phenotype (basename) and then Yes or No depending on whether the harmonisation was successful or not and in brackets the percentage of variants successfully harmonised
#You need to iterate over all the '.running.log' files to do this

rule summarise_harmonisation:
    input:
        log_files = expand(config['base_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/final/{input_tsv_gz_basename}_gwas_ssf.running.log', input_tsv_gz_basename=input_tsv_gz_basenames)
    output:
        summary_file = config['desired_output_folder']+'harmonisation_summary.tsv'
    run:
        summary_data = []
        for log_file in input.log_files:
            basename = os.path.basename(log_file).split('_gwas_ssf')[0]
            with open(log_file, 'r') as f:
                lines = f.readlines()
                success = any("Result\tSUCCESS_HARMONIZATION" in line for line in lines)
                if success:
                    percentage_line =  [line for line in lines if 'sites successfully harmonised' in line]
                    summary_data.append(f"{basename}\tYes\t({percentage_line})")
                else:
                    summary_data.append(f"{basename}\tNo\t(N/A)")
        
        with open(output.summary_file, 'w') as f:
            f.write("Phenotype\tHarmonised\tPercentage\n")
            for line in summary_data:
                f.write(line + "\n")

#Rule 4 = Run GCTA-COJO to output independent SNPs from the harmonised GWAS summary statistics but only for certain summary statistic files as determined in the config.yaml

# rule gcta_cojo:
#     input:
#         harmonised_gwas_summstats=lambda wildcards: config['desired_output_folder']+wildcards.input_tsv_gz_basename+'_gwas_ssf/final/'+wildcards.input_tsv_gz_basename+'_gwas_ssf.h.tsv.gz'
#     output:
#         independent_snps = config['desired_output_folder']+'{input_tsv_gz_basename}_gwas_ssf/final/{input_tsv_gz_basename}.jma'
#     resources:
#         mem_mb = 8000,
#         threads = 8
#     shell:
#         '''
#         1_gcta_cojo.sh {input.harmonised_gwas_summstats} {output.independent_snps}
#         '''