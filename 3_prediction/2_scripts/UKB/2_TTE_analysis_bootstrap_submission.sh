#!/bin/bash

# --- Define variable parameters here ---
# PREDICTORS="sex_Male,bmi,dbp,HF_prevalent"
# MODEL_NAME="MRpreds"
# REGULARISATION="" # Options: "ridge", "lasso", "ElasticNet", or "" (for unregularised)

# PREDICTORS='sex_Male,eGFR,AF_prevalent,CAD_prevalent,stroke_prevalent,T2D_prevalent,smoking_initiation'
# MODEL_NAME="nonMRpreds"
# REGULARISATION="" # Options: "ridge", "lasso", "ElasticNet", or "" (for unregularised)

# PREDICTORS='sex_Male,bmi,dbp,HF_prevalent,eGFR,AF_prevalent,CAD_prevalent,stroke_prevalent,T2D_prevalent,smoking_initiation'
# MODEL_NAME="allpreds"
# REGULARISATION="" # Options: "ridge", "lasso", "ElasticNet", or "" (for unregularised)
# REGULARISATION="ridge"
# REGULARISATION="lasso"
# REGULARISATION="ElasticNet"

#------------------------------------------------------------------------------------------------------
# --- Define required parameters ---
PROJECT_ID=$(dx env --bash | grep DX_PROJECT_CONTEXT_ID | cut -d '=' -f 2) #Retrieve the project ID

INPUT_TSV="${PROJECT_ID}":/jchan/MR_Prediction/3_output/1_DataProcessing/ukb_survival_df.tsv
INPUT_TSV_NAME="ukb_survival_df.tsv"
PYTHON_SCRIPT="${PROJECT_ID}":/jchan/MR_Prediction/2_scripts/2_TTE_analysis_bootstrap.py
PYTHON_SCRIPT_NAME="2_TTE_analysis_bootstrap.py"

OUTPUT_FOLDER="/jchan/MR_Prediction/3_output/2_TTE_Analysis/bootstrap/"

# --- Specify your Docker Image from a public registry (e.g., Docker Hub) ---
# REPLACE this with the name of the image you pushed
DOCKER_IMAGE="${PROJECT_ID}":/common/py312_sklearn_sksurv.tar.gz

# --- Optional Parameters (will use defaults in Python script if not set) ---
N_BOOTSTRAPS=100 # Default is 1 (for test purposes)
K_FOLDS=5     # Default is 5

OUTPUT_TSV_NAME="bootstrap_${N_BOOTSTRAPS}bs_${K_FOLDS}fold_${REGULARISATION}_${MODEL_NAME}.tsv"

# --- Construct the command line arguments ---
CMD_ARGS="--filename ${INPUT_TSV_NAME} "
CMD_ARGS+="--preds ${PREDICTORS} "
CMD_ARGS+="--output_filename ${OUTPUT_TSV_NAME} "

# Add optional arguments only if they are set
if [ -n "$N_BOOTSTRAPS" ]; then
    CMD_ARGS+="--n_bootstraps ${N_BOOTSTRAPS} "
fi
if [ -n "$K_FOLDS" ]; then
    CMD_ARGS+="--k_folds ${K_FOLDS} "
fi
if [ -n "$REGULARISATION" ]; then
    CMD_ARGS+="--regularisation ${REGULARISATION}"
fi

#------------------------------------------------------------------------------------------------------
# --- Construct the dx run command with the Docker image ---
dx run swiss-army-knife \
    -iin="${INPUT_TSV}" \
    -iin="${PYTHON_SCRIPT}" \
    -iimage_file="${DOCKER_IMAGE}" \
    -icmd="python ${PYTHON_SCRIPT_NAME} ${CMD_ARGS}" --instance-type "mem1_ssd1_v2_x4" --brief -y --destination "${OUTPUT_FOLDER}"