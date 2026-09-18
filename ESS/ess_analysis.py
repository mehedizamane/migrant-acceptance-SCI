import os
import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf

print("Loading datasets...")

# 1. LOAD ESS ROUND 7 HISTORICAL BASELINE
ess_raw = pd.read_csv("/kaggle/input/datasets/smmehedizaman/ess-r7/ESS7e02_3.csv", low_memory=False)

# Harmonize country codes to match standard ISO keys
ess_raw['cntry'] = ess_raw['cntry'].replace({'UK': 'GB'})

iso_mapping = {
    'AT': 'AUT', 'BE': 'BEL', 'CH': 'CHE', 'CZ': 'CZE', 'DE': 'DEU', 
    'DK': 'DNK', 'EE': 'EST', 'ES': 'ESP', 'FI': 'FIN', 'FR': 'FRA', 
    'GB': 'GBR', 'HU': 'HUN', 'IE': 'IRL', 'LT': 'LTU', 'NL': 'NLD', 
    'NO': 'NOR', 'PL': 'POL', 'PT': 'PRT', 'SE': 'SWE', 'SI': 'SVN'
}

# Filter master set down strictly to the 20 target tracking cohort countries
model_df = ess_raw[ess_raw['cntry'].isin(iso_mapping.keys())].copy()

# 2. RECODE INDIVIDUAL LEVEL VARIABLES (ESS Standard Layout)
print("Recoding individual level metrics...")

# A. Migration Acceptance Index (Switched to Impact Scale: 0 to 10)
# Variables: Economy (imbgeco), Culture (imueclt), Place to live (imwbcnt)
impact_cols = ['imbgeco', 'imueclt', 'imwbcnt']
for col in impact_cols:
    model_df[col] = pd.to_numeric(model_df[col], errors='coerce')
    # Remove missing codes (77: Refusal, 88: Don't know, 99: No answer)
    model_df.loc[model_df[col] > 10, col] = np.nan

# Calculate index as the mean across the three items
model_df['mig_accept_index'] = model_df[impact_cols].mean(axis=1)

# B. Intergroup Contact (z_contact using dfegcon)
model_df['contact_raw'] = pd.to_numeric(model_df['dfegcon'], errors='coerce')
model_df.loc[model_df['contact_raw'] > 7, 'contact_raw'] = np.nan
model_df['z_contact'] = (model_df['contact_raw'] - model_df['contact_raw'].mean()) / model_df['contact_raw'].std()

# C. Basic Demographics
model_df['age_years'] = pd.to_numeric(model_df['agea'], errors='coerce')
model_df.loc[model_df['age_years'] > 110, 'age_years'] = np.nan

model_df['female'] = pd.to_numeric(model_df['gndr'], errors='coerce').map({1: 0, 2: 1})
model_df['foreign_born'] = pd.to_numeric(model_df['brncntr'], errors='coerce').map({1: 0, 2: 1})

# D. Categorical Harmonizations (Education & Urbanicity)
model_df['edu_val'] = pd.to_numeric(model_df['eisced'], errors='coerce')
def map_edu(val):
    if val in [1, 2]: return 'Elementary_or_less'
    if val in [3, 4, 5]: return 'Secondary'
    if val in [6, 7]: return 'Tertiary'
    return np.nan
model_df['edu_cat'] = model_df['edu_val'].apply(map_edu)

model_df['urb_val'] = pd.to_numeric(model_df['domicil'], errors='coerce')
urb_mapping = {
    1: 'Big_city', 2: 'Suburbs_or_outskirts', 3: 'Town_or_small_city', 
    4: 'Country_village', 5: 'Farm_or_countryside'
}
model_df['urbanicity'] = model_df['urb_val'].map(urb_mapping)

# 3. FETCH WORLD BANK MACRO CONTROLS
print("Fetching indicators using robust 3-letter ISO keys...")
country_str_3 = ";".join(iso_mapping.values())

try:
    gdp_url = f"http://api.worldbank.org/v2/country/{country_str_3}/indicator/NY.GDP.PCAP.PP.CD?per_page=100&format=json&date=2013:2015"
    gdp_json = requests.get(gdp_url).json()[1]
    pop_url = f"http://api.worldbank.org/v2/country/{country_str_3}/indicator/SP.POP.TOTL?per_page=100&format=json&date=2013:2015"
    pop_json = requests.get(pop_url).json()[1]

    df_gdp = pd.DataFrame([{ 'iso3': item['countryiso3code'], 'gdp_pc': item['value'], 'gdp_year': item['date'] } for item in gdp_json if item['value'] is not None])
    df_pop = pd.DataFrame([{ 'iso3': item['countryiso3code'], 'pop': item['value'], 'pop_year': item['date'] } for item in pop_json if item['value'] is not None])

    df_gdp = df_gdp.sort_values('gdp_year').groupby('iso3').last().reset_index()
    df_pop = df_pop.sort_values('pop_year').groupby('iso3').last().reset_index()

    df_macro_raw = pd.merge(df_gdp, df_pop, on='iso3')
    reverse_mapping = {v: k for k, v in iso_mapping.items()}
    df_macro_raw['cntry'] = df_macro_raw['iso3'].map(reverse_mapping)

    df_macro_raw['log_population'] = np.log(df_macro_raw['pop'])
    df_macro_raw['gdp_pc_thousands'] = df_macro_raw['gdp_pc'] / 1000

    model_df = model_df.merge(df_macro_raw[['cntry', 'gdp_pc_thousands', 'log_population']], on='cntry', how='left')
except Exception as e:
    print(f"World Bank API fallback defaults used: {e}")
    model_df['gdp_pc_thousands'] = 45.0
    model_df['log_population'] = 16.0

# 4. DIRECTLY MAP FLAT GLOBAL SCI DATA FROM THE 4TH COLUMN
print("Processing flat Global SCI dataset...")
sci_df = pd.read_csv("/kaggle/input/datasets/smmehedizaman/global-sci-2020/global_sci_2020_2026.csv")

# Extract exactly the target elements using iloc positions to dodge string typing errors
sci_clean = pd.DataFrame({
    'cntry': sci_df.iloc[:, 0].astype(str).str.strip().replace({'UK': 'GB'}),
    'country_sci_raw': pd.to_numeric(sci_df.iloc[:, 3], errors='coerce')
})

# Merge clean metric structure directly into master set
model_df = model_df.merge(sci_clean, on='cntry', how='left')
model_df['z_country_sci'] = (model_df['country_sci_raw'] - model_df['country_sci_raw'].mean()) / model_df['country_sci_raw'].std()

# 5. PIPELINE CLEANING & QUALITY CHECK
required_columns = [
    'mig_accept_index', 'z_country_sci', 'z_contact', 
    'age_years', 'female', 'foreign_born', 'edu_cat', 'urbanicity',
    'gdp_pc_thousands', 'log_population'
]

model_df_clean = model_df.dropna(subset=required_columns).copy()

print("\n--- Data Attrition Complete ---")
print(f"Final Analysis Sample Size: {len(model_df_clean)} rows across {model_df_clean['cntry'].nunique()} countries.")

# 6. RUN STRUCTURAL MIXED-EFFECTS MODELS
ind_controls = (
    "age_years + female + foreign_born + "
    "C(edu_cat, Treatment('Elementary_or_less')) + "
    "C(urbanicity, Treatment('Farm_or_countryside'))"
)

f_m1 = f"mig_accept_index ~ z_country_sci + {ind_controls}"
f_m2 = f"mig_accept_index ~ z_contact + {ind_controls}"
f_m3 = f"mig_accept_index ~ z_country_sci + z_contact + {ind_controls}"
f_m4 = f"mig_accept_index ~ z_country_sci * z_contact + {ind_controls}"
f_m5 = f"mig_accept_index ~ z_country_sci * z_contact + gdp_pc_thousands + log_population + {ind_controls}"

models = [f_m1, f_m2, f_m3, f_m4, f_m5]
results = []

for idx, formula in enumerate(models, 1):
    print(f"Fitting Model {idx}...")
    res = smf.mixedlm(formula, data=model_df_clean, groups=model_df_clean['cntry']).fit(reml=True)
    results.append(res)

# Print execution report summaries
for idx, res in enumerate(results, 1):
    print(f"\n=== MODEL {idx} SUMMARY ===")
    print(res.summary())
