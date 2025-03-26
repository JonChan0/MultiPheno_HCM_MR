''' Snakemake script to run UNIdirectional Mendelian Randomisation 
Author: Jonathan Chan
Date: 2025-02-20
'''

rule all:
    input: 
     config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_mainline.tsv',
     config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_sensitivity_oddsratio.tsv'

rule direction1_mr:
    input:
        exposure_path = config['d1_exposure_path'],
        outcome_path = config['d1_outcome_path']
    params:
        selected_instruments = config['d1_instruments'],
        output_path = config['output_path'],
        exposure = config['d1_exposure'],
        outcome = config['d1_outcome'],
        ld_clumping = config['d1_ld_clump'],
        ld_eur_bed_file = config['ld_eur_bed_file'],
        plink_binary_path = config['plink_binary_path']
    resources:
        mem_mb = 32000
    output:
        tsv_output = config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_mainline.tsv',
        tsv2_output = config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_sensitivity_oddsratio.tsv'
    shell:
        '''
        module purge 
        module use -a /apps/eb/2022b/skylake/modules/all
        module load R/4.2.2-foss-2022b
        Rscript 2_TwoSampleMR.R {input.exposure_path} {input.outcome_path} {params.selected_instruments} {params.output_path} {params.exposure} {params.outcome} {params.ld_clumping} {params.ld_eur_bed_file} {params.plink_binary_path}
        '''