import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

os.makedirs("results/portability_analysis", exist_ok=True)

# 1. Figure ML1: Pooled CV versus leave-one-cohort-out performance for each model.
df_pooled = pd.read_csv("results/portability_analysis/model_performance_pooled_cv.csv")

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
models = df_pooled['Model'].values
x = np.arange(len(models))
width = 0.35

ax[0].bar(x - width/2, df_pooled['Mean_AUROC'], width, label='Pooled CV')
ax[0].bar(x + width/2, df_pooled['Mean_LOCO_AUROC'], width, label='Mean LOCO')
ax[0].set_ylabel('AUROC')
ax[0].set_title('AUROC Comparison')
ax[0].set_xticks(x)
ax[0].set_xticklabels(models)
ax[0].legend()
ax[0].set_ylim(0.4, 1.0)

ax[1].bar(x - width/2, df_pooled['Mean_AUPRC'], width, label='Pooled CV')
ax[1].bar(x + width/2, df_pooled['Mean_LOCO_AUPRC'], width, label='Mean LOCO')
ax[1].set_ylabel('AUPRC')
ax[1].set_title('AUPRC Comparison')
ax[1].set_xticks(x)
ax[1].set_xticklabels(models)
ax[1].legend()

plt.tight_layout()
plt.savefig("results/portability_analysis/Figure_ML1_Pooled_vs_LOCO.png", dpi=300)
plt.savefig("results/portability_analysis/Figure_ML1_Pooled_vs_LOCO.pdf")
plt.close()

# 2. Figure ML2: Held-out cohort performance
df_loco = pd.read_csv("results/portability_analysis/model_performance_leave_one_cohort_out.csv")
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
sns.barplot(data=df_loco, x='Model', y='AUROC', hue='Held_Out_Cohort', ax=ax[0])
ax[0].set_title('LOCO AUROC by Cohort')
ax[0].set_ylim(0.4, 1.0)
sns.barplot(data=df_loco, x='Model', y='AUPRC', hue='Held_Out_Cohort', ax=ax[1])
ax[1].set_title('LOCO AUPRC by Cohort')
plt.tight_layout()
plt.savefig("results/portability_analysis/Figure_ML2_LOCO_by_Cohort.png", dpi=300)
plt.savefig("results/portability_analysis/Figure_ML2_LOCO_by_Cohort.pdf")
plt.close()

# 3. Figure ML3: Biomarker Portability Gap
df_bpg = pd.read_csv("results/portability_analysis/biomarker_portability_gap.csv")
df_bpg_melt = df_bpg.melt(id_vars=['Model'], var_name='Metric', value_name='Gap')
plt.figure(figsize=(10, 5))
sns.barplot(data=df_bpg_melt, x='Model', y='Gap', hue='Metric')
plt.title('Biomarker Portability Gap')
plt.ylabel('Gap (Pooled - LOCO)')
plt.tight_layout()
plt.savefig("results/portability_analysis/Figure_ML3_BPG.png", dpi=300)
plt.savefig("results/portability_analysis/Figure_ML3_BPG.pdf")
plt.close()

# 4. Figure ML4: PCA before and after correction
df_raw = pd.read_csv("results/portability_analysis/dataset_clr.csv")
df_res = pd.read_csv("results/portability_analysis/dataset_clr_global_residualized.csv")

feature_cols = [c for c in df_raw.columns if c not in ['SampleID', 'Cohort', 'PD']]

pca_raw = PCA(n_components=2).fit_transform(df_raw[feature_cols].values)
pca_res = PCA(n_components=2).fit_transform(df_res[feature_cols].values)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
sns.scatterplot(x=pca_raw[:,0], y=pca_raw[:,1], hue=df_raw['Cohort'], style=df_raw['PD'], ax=ax[0], alpha=0.7)
ax[0].set_title('PCA Raw CLR')
ax[0].set_xlabel('PC1')
ax[0].set_ylabel('PC2')

sns.scatterplot(x=pca_res[:,0], y=pca_res[:,1], hue=df_res['Cohort'], style=df_res['PD'], ax=ax[1], alpha=0.7)
ax[1].set_title('PCA Global Residualized CLR')
ax[1].set_xlabel('PC1')
ax[1].set_ylabel('PC2')
plt.tight_layout()
plt.savefig("results/portability_analysis/Figure_ML4_PCA.png", dpi=300)
plt.savefig("results/portability_analysis/Figure_ML4_PCA.pdf")
plt.close()

# 5. Figure ML5: Cohort-Disease Dominance Ratio and disease-vector angle summary
# Assuming CDDR values and angle values are in CSVs
# Let's create a text-based plot or simple bar plot
cddr_df = pd.read_csv("results/portability_analysis/cohort_disease_dominance_ratio.csv")
geom_df = pd.read_csv("results/portability_analysis/disease_vector_geometry.csv")

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].bar(['Cohort R2', 'PD R2'], cddr_df['R2'].values, color=['blue', 'red'])
ax[0].set_title(f"PERMANOVA R2 (CDDR: {cddr_df['CDDR'].values[0]:.1f})")
ax[0].set_ylabel('R2 Explained Variance')

metrics = geom_df['Metric'].values
vals = geom_df['Value'].values
ax[1].bar(metrics, vals, color='purple')
ax[1].set_title("Disease Vector Geometry")
plt.tight_layout()
plt.savefig("results/portability_analysis/Figure_ML5_CDDR_Geometry.png", dpi=300)
plt.savefig("results/portability_analysis/Figure_ML5_CDDR_Geometry.pdf")
plt.close()

# 6. Figure ML6: Feature importance stability across cohorts
try:
    df_stab = pd.read_csv("results/portability_analysis/feature_stability_summary.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    df_stab_melt = df_stab.melt(id_vars=['Model'], var_name='Metric', value_name='Value')
    sns.barplot(data=df_stab_melt, x='Model', y='Value', hue='Metric')
    plt.title('Feature Stability Across Cohorts')
    plt.tight_layout()
    plt.savefig("results/portability_analysis/Figure_ML6_Feature_Stability.png", dpi=300)
    plt.savefig("results/portability_analysis/Figure_ML6_Feature_Stability.pdf")
    plt.close()
except:
    pass

print("Figures generated successfully.")
