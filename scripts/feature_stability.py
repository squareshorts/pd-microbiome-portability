import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
import xgboost as xgb
from sklearn.metrics import jaccard_score

SEED = 42

df = pd.read_csv("results/portability_analysis/dataset_clr.csv")
cohorts = df['Cohort'].unique()
feature_cols = [c for c in df.columns if c not in ['SampleID', 'Cohort', 'PD']]
X = df[feature_cols].values
y = df['PD'].values
cohort_labels = df['Cohort'].values

models = {
    'ElasticNet': LogisticRegressionCV(cv=5, penalty='elasticnet', solver='saga', l1_ratios=[0.1, 0.5, 0.9], max_iter=1000, random_state=SEED, n_jobs=-1, class_weight='balanced'),
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=SEED, class_weight='balanced', n_jobs=-1),
    'XGBoost': xgb.XGBClassifier(n_estimators=100, use_label_encoder=False, eval_metric='logloss', random_state=SEED, scale_pos_weight=(len(y)-sum(y))/sum(y), n_jobs=-1)
}

importance_records = []

for model_name, model in models.items():
    print(f"Feature extraction for {model_name}...")
    
    top10_sets = []
    top20_sets = []
    
    for held_out in cohorts:
        train_mask = cohort_labels != held_out
        X_train, y_train = X[train_mask], y[train_mask]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        
        model.fit(X_train_scaled, y_train)
        
        if model_name == 'ElasticNet':
            importances = model.coef_[0]
            abs_importances = np.abs(importances)
        else:
            importances = model.feature_importances_
            abs_importances = importances
            
        top10_idx = np.argsort(abs_importances)[-10:]
        top20_idx = np.argsort(abs_importances)[-20:]
        
        top10_sets.append(set(top10_idx))
        top20_sets.append(set(top20_idx))
        
        for i, f_col in enumerate(feature_cols):
            importance_records.append({
                'Model': model_name,
                'Held_Out_Cohort': held_out,
                'Feature': f_col,
                'Importance': importances[i],
                'Abs_Importance': abs_importances[i],
                'Rank': len(feature_cols) - np.where(np.argsort(abs_importances) == i)[0][0]
            })

pd.DataFrame(importance_records).to_csv("results/portability_analysis/feature_importance_by_model_and_cohort.csv", index=False)

# Compute stability summary
stability_summary = []
for model_name in models.keys():
    records = [r for r in importance_records if r['Model'] == model_name]
    
    # We have 3 cohorts, so 3 LOCO sets. We compute pairwise overlaps.
    held_outs = cohorts
    
    overlaps_10 = []
    overlaps_20 = []
    jaccards_20 = []
    
    for i in range(len(held_outs)):
        for j in range(i+1, len(held_outs)):
            ho1 = held_outs[i]
            ho2 = held_outs[j]
            
            set1_10 = set([r['Feature'] for r in records if r['Held_Out_Cohort'] == ho1 and r['Rank'] <= 10])
            set2_10 = set([r['Feature'] for r in records if r['Held_Out_Cohort'] == ho2 and r['Rank'] <= 10])
            overlaps_10.append(len(set1_10.intersection(set2_10)))
            
            set1_20 = set([r['Feature'] for r in records if r['Held_Out_Cohort'] == ho1 and r['Rank'] <= 20])
            set2_20 = set([r['Feature'] for r in records if r['Held_Out_Cohort'] == ho2 and r['Rank'] <= 20])
            overlaps_20.append(len(set1_20.intersection(set2_20)))
            
            jaccard = len(set1_20.intersection(set2_20)) / len(set1_20.union(set2_20))
            jaccards_20.append(jaccard)
            
    stability_summary.append({
        'Model': model_name,
        'Mean_Overlap_Top10': np.mean(overlaps_10),
        'Mean_Overlap_Top20': np.mean(overlaps_20),
        'Mean_Jaccard_Top20': np.mean(jaccards_20)
    })

pd.DataFrame(stability_summary).to_csv("results/portability_analysis/feature_stability_summary.csv", index=False)
print("Feature stability analysis complete.")
