import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import statsmodels.api as sm

# Load CLR data
df = pd.read_csv("results/portability_analysis/dataset_clr.csv")

feature_cols = [c for c in df.columns if c not in ['SampleID', 'Cohort', 'PD']]
X = df[feature_cols].values

# 1. PC1 Loading vector
pca = PCA(n_components=1)
pca.fit(X)
pc1_loading = pca.components_[0] # length = num_features

# 2. Cohort-adjusted PD coefficient vector
# Formula: clr_feature ~ PD + Cohort
# We will fit an OLS for each feature

# Dummy code cohort (drop one to avoid collinearity)
cohort_dummies = pd.get_dummies(df['Cohort'], drop_first=True)
exog = pd.concat([df[['PD']], cohort_dummies], axis=1)
exog = sm.add_constant(exog)

pd_effects = []
for col in feature_cols:
    endog = df[col]
    model = sm.OLS(endog, exog.astype(float)).fit()
    pd_effects.append(model.params['PD'])

pd_effects_vector = np.array(pd_effects)

# 3. Geometry
def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_sim = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_sim))
    return cos_sim, angle

cos_sim, angle = angle_between(pc1_loading, pd_effects_vector)
# If angle > 90, sometimes we care about the absolute alignment (so we might check 180 - angle)
# But let's just report the direct angle

print(f"Cosine similarity: {cos_sim:.3f}")
print(f"Angle: {angle:.1f} degrees")
print(f"PC1 variance explained: {pca.explained_variance_ratio_[0]*100:.1f}%")

res_df = pd.DataFrame([{
    "Metric": "Cosine_Similarity", "Value": cos_sim
}, {
    "Metric": "Angle_Degrees", "Value": angle
}, {
    "Metric": "PC1_Variance_Explained_Pct", "Value": pca.explained_variance_ratio_[0]*100
}])

res_df.to_csv("results/portability_analysis/disease_vector_geometry.csv", index=False)
