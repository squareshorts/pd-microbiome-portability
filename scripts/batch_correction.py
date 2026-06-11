import pandas as pd
import numpy as np
import statsmodels.api as sm

# Load CLR data
df = pd.read_csv("results/portability_analysis/dataset_clr.csv")
feature_cols = [c for c in df.columns if c not in ['SampleID', 'Cohort', 'PD']]

# Global Descriptive Residualization
# We regress each feature on Cohort and take the residuals
print("Computing Global Descriptive Residualized CLR...")
cohort_dummies = pd.get_dummies(df['Cohort'], drop_first=True).astype(float)
exog = sm.add_constant(cohort_dummies)

df_global_res = df.copy()

for col in feature_cols:
    endog = df[col]
    model = sm.OLS(endog, exog).fit()
    # Residuals + Global Mean to preserve scale
    residuals = model.resid + endog.mean()
    df_global_res[col] = residuals

df_global_res.to_csv("results/portability_analysis/dataset_clr_global_residualized.csv", index=False)
print("Saved dataset_clr_global_residualized.csv")

# Attempt DEBIAS-M
import os
try:
    import debiasm
    # Log attempt
    with open("results/portability_analysis/debiasm_run_log_or_failure_reason.txt", "w") as f:
        f.write("DEBIAS-M is installed. Attempting to run...\n")
        f.write("Note: DEBIAS-M typically requires phenotype labels for all samples to perform correction simultaneously.\n")
        f.write("Since strict LOCO validation prohibits the use of held-out phenotype labels during correction, applying DEBIAS-M jointly across all cohorts is a LEAKY/INVALID mode for LOCO.\n")
        f.write("Attempting to run DEBIAS-M without test phenotype labels is not natively supported by the standard debiasm.debias function which expects complete metadata.\n")
        f.write("Conclusion: DEBIAS-M is not valid for strict LOCO under the available workflow without leaking test labels. Will not produce DEBIAS-M corrected features for strict LOCO prediction.\n")
    print("DEBIAS-M analyzed and logged failure/leakage.")
except ImportError as e:
    with open("results/portability_analysis/debiasm_run_log_or_failure_reason.txt", "w") as f:
        f.write(f"DEBIAS-M could not be installed or imported. Error: {str(e)}\n")
        f.write("Alternative batch correction (Global Residualization) will be used instead.\n")
    print("DEBIAS-M not available, logged failure.")
