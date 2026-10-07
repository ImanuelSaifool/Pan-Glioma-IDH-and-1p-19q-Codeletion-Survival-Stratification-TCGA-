# Imports
import os

import pandas as pd
import matplotlib.pyplot as plt
import os
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.statistics import pairwise_logrank_test
from lifelines import CoxPHFitter
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# lgg reading and cleaning
Igg_df = pd.read_csv("C:\\Users\\imanu\\Downloads\\lgg_tcga_pan_can_atlas_2018\\lgg_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", sep='\t', comment='#')
Igg_df = Igg_df[Igg_df['IN_PANCANPATHWAYS_FREEZE'] == 'Yes']
Igg_df.drop(columns=['FORM_COMPLETION_DATE','INFORMED_CONSENT_VERIFIED','OTHER_PATIENT_ID','ICD_O_3_HISTOLOGY', 'PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT', 'PATH_N_STAGE', 'PATH_M_STAGE', 'PATH_T_STAGE', 'AJCC_STAGING_EDITION', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'WEIGHT', 'DFS_MONTHS', 'DFS_STATUS', 'NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT', 'PRIOR_DX'], inplace=True)
Igg_df['OS_STATUS'] = Igg_df['OS_STATUS'].str[0].astype('Int64')
Igg_df['DSS_STATUS'] = Igg_df['DSS_STATUS'].str[0].astype('Int64')
Igg_df['PFS_STATUS'] = Igg_df['PFS_STATUS'].str[0].astype('Int64')
Igg_df['IDH_STATUS'] = Igg_df['SUBTYPE'].str.extract(r'(IDHwt|IDHmut-codel|IDHmut-non-codel)')
Igg_df = Igg_df.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# gbm reading and cleaning
Gdf = pd.read_csv("C:\\Users\\imanu\\Downloads\\gbm_tcga_pan_can_atlas_2018\\gbm_tcga_pan_can_atlas_2018\\data_clinical_patient.txt", sep='\t', comment='#')
Gdf = Gdf[Gdf['IN_PANCANPATHWAYS_FREEZE'] == 'Yes']
Gdf.drop(columns=['FORM_COMPLETION_DATE','INFORMED_CONSENT_VERIFIED','OTHER_PATIENT_ID','ICD_O_3_HISTOLOGY', 'PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT', 'PATH_N_STAGE', 'PATH_M_STAGE', 'PATH_T_STAGE', 'AJCC_STAGING_EDITION', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'WEIGHT', 'DFS_MONTHS', 'DFS_STATUS', 'NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT', 'PRIOR_DX'], inplace=True)
Gdf['IDH_STATUS'] = Gdf['SUBTYPE'].str.extract(r'(IDHwt|IDHmut-non-codel)')
Gdf.loc[Gdf['SUBTYPE'] == 'GBM', 'IDH_STATUS'] = 'IDHwt'
cols = ['AGE', 'SEX', 'DAYS_TO_BIRTH', 'ICD_O_3_SITE', 'ICD_10', 
        'HISTORY_NEOADJUVANT_TRTYN', 'DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS']
Gdf = Gdf.drop(Gdf[Gdf[cols].isna().all(axis=1)].index)
Gdf['OS_STATUS'] = Gdf['OS_STATUS'].str[0].astype('Int64')
Gdf['PFS_STATUS'] = Gdf['PFS_STATUS'].str[0].astype('Int64')
Gdf['DSS_STATUS'] = Gdf['DSS_STATUS'].str[0].astype('Int64')
Gdf = Gdf.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combined data file
combined_df = pd.concat([Igg_df, Gdf], ignore_index=True)
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Figures
kmf = KaplanMeierFitter()
kmf.fit(durations=combined_df['OS_MONTHS'], event_observed=combined_df['OS_STATUS'], label='All Patients')
kmf.plot_survival_function()
plt.xlabel('Months')
plt.ylabel('Survival probability')
plt.title('Overall survival: all patients')
plt.show()

os.makedirs("figures", exist_ok=True)
fig, ax = plt.subplots(figsize=(10, 6))
for g, d in combined_df.dropna(subset=['IDH_STATUS']).groupby('IDH_STATUS'):
    KaplanMeierFitter().fit(d['OS_MONTHS'], d['OS_STATUS'], label=g).plot_survival_function(ax=ax)
ax.set_xlabel('Months')
ax.set_ylabel('Overall survival probability')
ax.set_title('Overall survival by IDH/1p19q subtype (TCGA)')
fig.savefig("figures/km_by_subtype.png", dpi=300, bbox_inches="tight")
plt.show()
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results = pairwise_logrank_test(
    event_durations=combined_df['OS_MONTHS'],
    groups=combined_df['IDH_STATUS'],
    event_observed=combined_df['OS_STATUS']
)
print(results.summary.to_string())
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
features = ['OS_MONTHS', 'OS_STATUS', 'AGE', 'IDH_STATUS']
cox_df = combined_df[features].copy()

# 2. Convert text categories (like SUBTYPE) into numeric 1/0 dummy variables
cox_df = pd.get_dummies(cox_df, columns=['IDH_STATUS'], drop_first=True)

# 3. Drop rows with missing values in these specific columns
cox_df = cox_df.dropna()

# 4. Initialize and fit the Cox model using ONLY this clean subset
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='OS_MONTHS', event_col='OS_STATUS')
cph.print_summary()

for g, d in combined_df.dropna(subset=['IDH_STATUS']).groupby('IDH_STATUS'):
    k = KaplanMeierFitter().fit(d['OS_MONTHS'], d['OS_STATUS'])
    print(g, len(d), int(d['OS_STATUS'].sum()), round(k.median_survival_time_, 1))
results.print_summary()

cph.check_assumptions(cox_df, p_value_threshold=0.05, show_plots=True)

