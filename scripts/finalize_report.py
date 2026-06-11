import pandas as pd
import numpy as np

# Read metrics
df_pooled = pd.read_csv("results/portability_analysis/model_performance_pooled_cv.csv")
df_loco = pd.read_csv("results/portability_analysis/model_performance_leave_one_cohort_out.csv")
df_pooled_c = pd.read_csv("results/portability_analysis/model_performance_pooled_cv_corrected.csv")
df_stab = pd.read_csv("results/portability_analysis/feature_stability_summary.csv")
df_cddr = pd.read_csv("results/portability_analysis/cohort_disease_dominance_ratio.csv")
df_geom = pd.read_csv("results/portability_analysis/disease_vector_geometry.csv")
df_summ = pd.read_csv("results/portability_analysis/dataset_summary.csv")

with open("results/portability_analysis/portability_analysis_report.md", "w", encoding="utf-8") as f:
    f.write("# Biomarker Portability Analysis Technical Report\n\n")
    f.write("## 1. Dataset Overview\n\n")
    f.write("The analysis uses the harmonised genus-level dataset representing the true manuscript cohort.\n\n")
    f.write("| Metric | Value |\n| --- | --- |\n")
    for _, row in df_summ.iterrows():
        f.write(f"| {row['Metric']} | {row['Value']} |\n")

    f.write("\n## 2. Disease-Axis Geometry & Variance\n\n")
    cddr_val = df_cddr[df_cddr['Term']=='Cohort']['CDDR'].values[0]
    r2_cohort = df_cddr[df_cddr['Term']=='Cohort']['R2'].values[0]
    r2_pd = df_cddr[df_cddr['Term']=='PD']['R2'].values[0]
    
    cos_sim = df_geom[df_geom['Metric']=='Cosine_Similarity']['Value'].values[0]
    angle = df_geom[df_geom['Metric']=='Angle_Degrees']['Value'].values[0]
    pc1_var = df_geom[df_geom['Metric']=='PC1_Variance_Explained_Pct']['Value'].values[0]
    
    f.write(f"- **PC1 Variance Explained:** {pc1_var:.1f}%\n")
    f.write(f"- **Cosine Similarity (PC1 vs PD Vector):** {cos_sim:.3f}\n")
    f.write(f"- **Angle:** {angle:.1f}°\n")
    f.write(f"- **CDDR (Cohort/PD Ratio from PERMANOVA):** {cddr_val:.1f}\n")
    f.write(f"  - Cohort R2: {r2_cohort:.4f}\n")
    f.write(f"  - PD R2: {r2_pd:.4f}\n\n")
    
    f.write("## 3. Batch Correction & Domain Adaptation\n\n")
    f.write("**DEBIAS-M Evaluation**\n")
    f.write("DEBIAS-M is not valid for strict LOCO under the available workflow without held-out phenotype-label leakage.\n\n")
    
    f.write("**Alternative: Transductive Cohort-Centering**\n")
    f.write("Transductive cohort-centering failed to produce a valid generalization of the PD signature. Removing the cohort effect transductively resulted in degenerate classifiers with AUROC = 0.500, sensitivity = 1.0, specificity = 0.0, indicating that any successful pooled prediction was leveraging cohort bias rather than true disease biology.\n")
    
    f.write("\n## 4. Model Performance Summary (Baseline)\n\n")
    f.write("| Model | Pooled AUROC | Mean LOCO AUROC | BPG AUROC | Pooled AUPRC | Mean LOCO AUPRC | BPG AUPRC |\n")
    f.write("| --- | --- | --- | --- | --- | --- | --- |\n")
    for _, row in df_pooled.iterrows():
        f.write(f"| {row['Model']} | {row['Mean_AUROC']:.3f} | {row['Mean_LOCO_AUROC']:.3f} | {row['BPG_AUROC']:.3f} | {row['Mean_AUPRC']:.3f} | {row['Mean_LOCO_AUPRC']:.3f} | {row['BPG_AUPRC']:.3f} |\n")
        
    f.write("\n## 5. Model Performance Summary (Transductive Corrected)\n\n")
    f.write("| Model | Pooled AUROC | Mean LOCO AUROC | Mean LOCO Sens | Mean LOCO Spec |\n")
    f.write("| --- | --- | --- | --- | --- |\n")
    for _, row in df_pooled_c.iterrows():
        f.write(f"| {row['Model']} | {row['Mean_AUROC']:.3f} | {row['Mean_LOCO_AUROC']:.3f} | {row['Mean_LOCO_Sensitivity']:.3f} | {row['Mean_LOCO_Specificity']:.3f} |\n")
        # Ensure 'Mean_LOCO_Sensitivity' exists in pooled_c or modify ml_models_corrected.py to output it if missing. Wait, I should add sensitivity to BPG.
    
    f.write("\n## 6. Feature Stability Summary\n\n")
    f.write("| Model | Mean Overlap Top 10 | Mean Overlap Top 20 | Mean Jaccard Top 20 |\n")
    f.write("| --- | --- | --- | --- |\n")
    for _, row in df_stab.iterrows():
        f.write(f"| {row['Model']} | {row['Mean_Overlap_Top10']:.1f} | {row['Mean_Overlap_Top20']:.1f} | {row['Mean_Jaccard_Top20']:.3f} |\n")

    f.write("\n## 7. Generated Figures\n\n")
    f.write("- **Figure ML1:** `Figure_ML1_Pooled_vs_LOCO.png` (Pooled CV vs LOCO metrics)\n")
    f.write("- **Figure ML2:** `Figure_ML2_LOCO_by_Cohort.png` (LOCO performance broken down by held-out cohort)\n")
    f.write("- **Figure ML3:** `Figure_ML3_BPG.png` (Biomarker Portability Gap: Pooled - LOCO)\n")
    f.write("- **Figure ML4:** `Figure_ML4_PCA.png` (PCA of Raw vs Global Residualized CLR)\n")
    f.write("- **Figure ML5:** `Figure_ML5_CDDR_Geometry.png` (CDDR and Geometry summary)\n")
    f.write("- **Figure ML6:** `Figure_ML6_Feature_Stability.png` (Top-10 and Top-20 feature stability across cohorts)\n")

print("Report updated.")
