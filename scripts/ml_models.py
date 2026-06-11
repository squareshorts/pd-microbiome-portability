import os
import pandas as pd
import numpy as np
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, confusion_matrix, precision_recall_curve, auc
from sklearn.pipeline import Pipeline
import xgboost as xgb
import warnings
warnings.filterwarnings('ignore')

SEED = 42

def compute_metrics(y_true, y_pred_prob, threshold=0.5):
    y_pred = (y_pred_prob >= threshold).astype(int)
    auroc = roc_auc_score(y_true, y_pred_prob)
    auprc = average_precision_score(y_true, y_pred_prob)
    bac = balanced_accuracy_score(y_true, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    return auroc, auprc, bac, sensitivity, specificity

def get_optimal_threshold(y_true, y_pred_prob):
    # Select threshold that maximizes balanced accuracy on training set
    thresholds = np.linspace(0.01, 0.99, 99)
    best_thresh = 0.5
    best_bac = 0
    for t in thresholds:
        y_pred = (y_pred_prob >= t).astype(int)
        bac = balanced_accuracy_score(y_true, y_pred)
        if bac > best_bac:
            best_bac = bac
            best_thresh = t
    return best_thresh

# Load CLR data
df = pd.read_csv("results/portability_analysis/dataset_clr.csv")
cohorts = df['Cohort'].unique()

# Exclude metadata columns to get features
feature_cols = [c for c in df.columns if c not in ['SampleID', 'Cohort', 'PD']]
X = df[feature_cols].values
y = df['PD'].values
cohort_labels = df['Cohort'].values

# Models
models = {
    'ElasticNet': LogisticRegressionCV(cv=5, penalty='elasticnet', solver='saga', l1_ratios=[0.1, 0.5, 0.9], max_iter=1000, random_state=SEED, n_jobs=4, class_weight='balanced'),
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=SEED, class_weight='balanced', n_jobs=4),
    'XGBoost': xgb.XGBClassifier(n_estimators=100, use_label_encoder=False, eval_metric='logloss', random_state=SEED, scale_pos_weight=(len(y)-sum(y))/sum(y), n_jobs=4)
}

pooled_results = []
loco_results = []

print("Running ML evaluations...")

for model_name, model in models.items():
    print(f"--- {model_name} ---")
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', model)
    ])
    
    # 1. Pooled CV
    print("Running Pooled CV...")
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=5, random_state=SEED)
    
    aurocs, auprcs, bacs, sens, specs = [], [], [], [], []
    prevalences = []
    
    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        pipeline.fit(X_train, y_train)
        y_pred_prob = pipeline.predict_proba(X_test)[:, 1]
        
        # Determine threshold on train
        y_train_pred_prob = pipeline.predict_proba(X_train)[:, 1]
        thresh = get_optimal_threshold(y_train, y_train_pred_prob)
        
        auroc, auprc, bac, sensitivity, specificity = compute_metrics(y_test, y_pred_prob, thresh)
        
        aurocs.append(auroc)
        auprcs.append(auprc)
        bacs.append(bac)
        sens.append(sensitivity)
        specs.append(specificity)
        prevalences.append(np.mean(y_test))
        
    pooled_results.append({
        'Model': model_name,
        'Mean_AUROC': np.mean(aurocs),
        'Mean_AUPRC': np.mean(auprcs),
        'Mean_Balanced_Accuracy': np.mean(bacs),
        'Mean_Sensitivity': np.mean(sens),
        'Mean_Specificity': np.mean(specs),
        'Mean_Prevalence': np.mean(prevalences)
    })
    
    # 2. Leave-One-Cohort-Out
    print("Running LOCO CV...")
    loco_aurocs, loco_auprcs, loco_bacs = [], [], []
    for held_out in cohorts:
        train_mask = cohort_labels != held_out
        test_mask = cohort_labels == held_out
        
        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]
        
        pipeline.fit(X_train, y_train)
        y_pred_prob = pipeline.predict_proba(X_test)[:, 1]
        
        # Determine threshold on train
        y_train_pred_prob = pipeline.predict_proba(X_train)[:, 1]
        thresh = get_optimal_threshold(y_train, y_train_pred_prob)
        
        auroc, auprc, bac, sensitivity, specificity = compute_metrics(y_test, y_pred_prob, thresh)
        
        loco_results.append({
            'Model': model_name,
            'Held_Out_Cohort': held_out,
            'AUROC': auroc,
            'AUPRC': auprc,
            'Balanced_Accuracy': bac,
            'Sensitivity': sensitivity,
            'Specificity': specificity,
            'Test_Prevalence': np.mean(y_test),
            'Test_N': len(y_test),
            'Test_Positives': sum(y_test)
        })
        loco_aurocs.append(auroc)
        loco_auprcs.append(auprc)
        loco_bacs.append(bac)
        
    # Calculate Mean LOCO for BPG
    mean_loco_auroc = np.mean(loco_aurocs)
    mean_loco_auprc = np.mean(loco_auprcs)
    mean_loco_bac = np.mean(loco_bacs)
    
    pooled_results[-1]['Mean_LOCO_AUROC'] = mean_loco_auroc
    pooled_results[-1]['Mean_LOCO_AUPRC'] = mean_loco_auprc
    pooled_results[-1]['Mean_LOCO_Balanced_Accuracy'] = mean_loco_bac
    pooled_results[-1]['BPG_AUROC'] = np.mean(aurocs) - mean_loco_auroc
    pooled_results[-1]['BPG_AUPRC'] = np.mean(auprcs) - mean_loco_auprc
    pooled_results[-1]['BPG_Balanced_Accuracy'] = np.mean(bacs) - mean_loco_bac

pd.DataFrame(pooled_results).to_csv("results/portability_analysis/model_performance_pooled_cv.csv", index=False)
pd.DataFrame(loco_results).to_csv("results/portability_analysis/model_performance_leave_one_cohort_out.csv", index=False)

bpg_df = pd.DataFrame(pooled_results)[['Model', 'BPG_AUROC', 'BPG_AUPRC', 'BPG_Balanced_Accuracy']]
bpg_df.to_csv("results/portability_analysis/biomarker_portability_gap.csv", index=False)

print("ML evaluations complete.")
