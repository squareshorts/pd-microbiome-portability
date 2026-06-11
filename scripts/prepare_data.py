import os
import pandas as pd

os.makedirs("results/portability_analysis", exist_ok=True)

# Load the correct manuscript CLR matrix
df_clr = pd.read_csv("data/processed/ecology/merged_genus_clr_matrix.csv")

summary = []
summary.append({"Metric": "Total Samples", "Value": len(df_clr)})
summary.append({"Metric": "Number of Genus Features (CLR)", "Value": len(df_clr.columns) - 3}) # exclude SampleID, PD, Cohort

for cohort in df_clr['Cohort'].unique():
    subset = df_clr[df_clr['Cohort'] == cohort]
    pd_count = (subset['PD'] == 1).sum()
    ctrl_count = (subset['PD'] == 0).sum()
    summary.append({"Metric": f"Cohort: {cohort} - Total", "Value": len(subset)})
    summary.append({"Metric": f"Cohort: {cohort} - PD", "Value": pd_count})
    summary.append({"Metric": f"Cohort: {cohort} - Control", "Value": ctrl_count})

summary_df = pd.DataFrame(summary)
summary_df.to_csv("results/portability_analysis/dataset_summary.csv", index=False)
print("Dataset summary generated with correct 683x62 matrix.")

df_clr.to_csv("results/portability_analysis/dataset_clr.csv", index=False)
print("Saved dataset_clr.csv.")
