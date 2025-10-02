import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, RandomizedSearchCV
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc

def fit_and_score_bootstrap_cv(filename, n_bootstraps, k_folds,preds, output_filename, regularisation,  ):
    """
    Performs bootstrap sampling and k-fold cross-validation on survival data
    with selectable regularisation methods.

    For each bootstrap sample, it performs k-fold cross-validation and computes
    various survival metrics. The results are saved to a TSV file.
    """
    try:
        df = pd.read_csv(filename, sep='\t')
    except FileNotFoundError:
        print(f"Error: Input file not found at {filename}")
        return

    predictors = preds.split(',')
    df.dropna(axis=0, subset=predictors, inplace=True)

    results = {}
    # Define the hyperparameter search space for Ridge
    ridge_alpha_params = {'alpha': np.logspace(-4, 2, 100)}
    all_eval_times = np.array([40, 50, 60, 70, 80, 90])

    for i in range(n_bootstraps):
        # Create bootstrap sample
        sample_df = df.sample(n=len(df), replace=True, random_state=i)

        # Prepare data for scikit-survival
        data_y = np.array(
            list(zip(sample_df['indicator'], sample_df['survival_age'])),
            dtype=[('indicator', bool), ('time', '<f8')]
        )

        data_x = sample_df[predictors].to_numpy()

        cv = KFold(n_splits=k_folds, shuffle=True, random_state=26)

        # *** FIX: Initialize all keys for fold_scores at the start ***
        fold_scores = {"Harrells_C": [], "Unos_C": []}
        for t in all_eval_times:
            fold_scores[f"AUC_{t}"] = []

        for train_idx, val_idx in cv.split(data_x):
            x_train, x_val = data_x[train_idx], data_x[val_idx]
            y_train, y_val = data_y[train_idx], data_y[val_idx]
            
            estimator = None # Ensure estimator is defined

            # --- Select, fit, and tune the model based on regularisation arg ---
            if regularisation is None:
                estimator = CoxPHSurvivalAnalysis()
                estimator.fit(x_train, y_train)

            elif regularisation == 'ridge':
                # Use RandomizedSearchCV for hyperparameter tuning within the training fold
                base_estimator = CoxPHSurvivalAnalysis()
                ridge_search = RandomizedSearchCV(
                    base_estimator,
                    param_distributions=ridge_alpha_params,
                    cv=3, # Inner CV loop for tuning
                    random_state=26,
                    n_iter=50, # Number of parameter settings to sample
                    n_jobs=-1
                )
                ridge_search.fit(x_train, y_train)
                estimator = ridge_search.best_estimator_

            elif regularisation == 'lasso':
                estimator = CoxnetSurvivalAnalysis(l1_ratio=1.0, alpha_min_ratio='auto', max_iter=100000)
                estimator.fit(x_train, y_train)

            elif regularisation == 'ElasticNet':
                estimator = CoxnetSurvivalAnalysis(l1_ratio=0.5, alpha_min_ratio='auto', max_iter=100000)
                estimator.fit(x_train, y_train)

            else:
                raise ValueError(f"Unknown regularisation type: {regularisation}")

            # --- Calculate Evaluation Metrics on the validation fold ---
            risk_scores = estimator.predict(x_val)
            
            # Harrell's C-index
            harrells_c = concordance_index_censored(y_val['indicator'], y_val['time'], risk_scores)[0]
            fold_scores["Harrells_C"].append(harrells_c)

            # Uno's C-index
            unos_c = concordance_index_ipcw(y_train, y_val, risk_scores)[0]
            fold_scores["Unos_C"].append(unos_c)

             # --- Time-dependent AUC calculation for individual time points ---
            min_val_time, max_val_time = y_val['time'].min(), y_val['time'].max()
            valid_eval_times = all_eval_times[(all_eval_times > min_val_time) & (all_eval_times < max_val_time)]
            
            auc_at_valid_times = {}
            if len(valid_eval_times) > 0:
                auc_values, _ = cumulative_dynamic_auc(y_train, y_val, risk_scores, valid_eval_times)
                for idx, t in enumerate(valid_eval_times):
                    auc_at_valid_times[t] = auc_values[idx]
            
            # Append scores (or NaN if not calculable) for all desired time points
            for t in all_eval_times:
                fold_scores[f"AUC_{t}"].append(auc_at_valid_times.get(t, np.nan))

        # --- Aggregate scores for the bootstrap iteration ---
        results[i] = {
            "Harrells_C": np.mean(fold_scores["Harrells_C"]),
            "Unos_C": np.mean(fold_scores["Unos_C"]),
        }
        for t in all_eval_times:
            results[i][f"AUC_{t}"] = np.nanmean(fold_scores[f"AUC_{t}"])
            
        print(f"Bootstrap {i+1}/{n_bootstraps} complete.")

    # Save results to a TSV file
    output_df = pd.DataFrame.from_dict(results, orient='index')
    output_df.index.name = 'bootstrap_iteration'
    output_df.to_csv(output_filename, sep='\t')
    print(f"\nSuccessfully saved results to {output_filename}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform bootstrap sampling and k-fold cross-validation for survival analysis.")
    parser.add_argument('--filename', type=str, required=True, help='Input TSV file path.')
    parser.add_argument('--n_bootstraps', type=int, default=1, help='Number of bootstrap samples.')
    parser.add_argument('--k_folds', type=int, default=5, help='The "k" for k-fold cross-validation.')
    parser.add_argument('--preds', required=True, type=str,help='List of predictors to use in the model.')
    parser.add_argument('--output_filename', type=str, required=True, help='Name for the output TSV file.')
    parser.add_argument('--regularisation', type=str, default=None, help="Type of regularisation. Options: 'ridge', 'lasso', 'ElasticNet'. Leave blank for unregularised.")

    args = parser.parse_args()
    
    # Handle the case where no regularisation is passed
    regularisation_arg = args.regularisation if args.regularisation and args.regularisation.lower() != 'none' else None


    fit_and_score_bootstrap_cv(
        args.filename,
        args.n_bootstraps,
        args.k_folds,
        args.preds,
        args.output_filename,
        regularisation_arg
    )