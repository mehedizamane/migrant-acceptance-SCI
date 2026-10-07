# Cross-border Social Connectedness and Migrant Acceptance Across Geographic Scales

Code for the Gallup World Poll and Meta Social Connectedness Index (SCI) analyses in
*Cross-border Social Connectedness and Migrant Acceptance Across Geographic Scales*
(Sehgal, Ohamadike, Lin, Zaman, Coimbra Vieira, Evans, Owoo, and Zagheni).

This directory reproduces **Figures 1, 2, 4, and 5** and **SI Tables 1 and 8–17**. SI Tables 2–3
describe variable coding, which is implemented in stages 02 and 04.
The European Social Survey and Eurobarometer analyses (Figure 3; SI Tables 4–7 and
18–27) were conducted separately by co-authors and are not included here.

Seven small public-source inputs are included under `data/raw/`; everything else must be
obtained separately. See [Inputs](#inputs).

## Pipeline

From this `gallup/` directory, run every stage in order with `python run_all.py`, or run any stage on its own.
Stages communicate only through files. The two R stages take their input and output
paths as arguments; `run_all.py` shows the exact commands.

| Stage | Produces | Paper |
|---|---|---|
| `00_build_gadm1_subset.py` | GADM 4.1 ADM1 layer from the full GADM download | |
| `01_build_sci_measures.py` | Country GlobalSCI (2020 and 2026 releases), ADM1 RegionalSCI, between/within variance decomposition, release correlation | Eq. 1, 3–4, 7; Fig. 2A–B |
| `02_build_controls.py` | National, ADM1, and CEPII dyadic controls | SI Table 3 |
| `03_build_geography_adjusted_sci.py` | Geography- and history-adjusted country and regional SCI | Eq. 5–6 |
| `04_build_gallup_analysis_data.py` | Four-wave Gallup analytic file: outcomes, controls, weights | SI Table 2 |
| `05_prepare_country_model_input.py` | Standardized model input and balanced country samples | |
| `06_fit_country_models.R` | Country-level mixed and OLS models | Fig. 1, 2D, 4, 5; SI Tables 8–12 (top panels), 16–17 |
| `07_seed_region_compositions.py` | Documented multi-region compositions | |
| `08_match_gallup_regions.py` | Candidate matches of Gallup region codes to GADM1 | |
| `09_link_gallup_regions_to_adm1.py` | Deterministic one-to-one ADM1 linkage; regional model input | SI Tables 1, 15 |
| `10_fit_regional_models.R` | Regional models | Fig. 2C; SI Tables 12 (bottom panel), 13–14 |
| `11_make_figures_and_tables.py` | Figures and one CSV per SI table | All of the above |

Outputs are written to `outputs/`: figures to `outputs/figures/`, and each SI table plus
the values quoted only in the text (`in_text_values.csv`) to `outputs/tables/`.

All regional analyses (Fig. 2C; SI Tables 12–15) use the deterministic one-to-one
Gallup-to-ADM1 linkage built by stages 08–09.

## Running

```bash
cd gallup
python run_all.py
```

The results in the paper were produced with the versions below. Other recent versions
should work, but have not been checked.

- **Python 3.9.21**: pandas 2.2.3, numpy 1.26.4, scipy 1.13.1, matplotlib 3.9.2,
  geopandas 1.0.1, pyogrio 0.10.0, pyproj 3.6.1, pyarrow 21.0.0, openpyxl 3.1.5,
  pycountry 24.6.1, unidecode 1.4.0
- **R 4.4.1** (`Rscript` on the PATH): data.table 1.17.0, lme4 1.1.37, lmerTest 3.1.3,
  fixest 0.12.1, Matrix 1.7.0
- The system `unzip` command, which stages 04, 08, and 09 use to stream the Gallup archive

A full run takes about 30 minutes (plus about 10 minutes for stage 00 on the first run)
and needs about 2.5 GB of RAM, 5.2 GB for the inputs, and about 1 GB for intermediates.
Stages 04, 08, and 09 stream the Gallup archive (about 17.5 GB uncompressed) without
extracting it; stage 08 caches region-code counts on its first run.

## Inputs

Place every input under `gallup/data/`; paths below are relative to it. Full SHA-256 checksums identify the exact files used.

### Included in this directory

| Path | Contents |
|---|---|
| `raw/controls/country_metadata_reference.csv` | ISO codes, UN numeric codes, country names, World Bank regions |
| `raw/controls/un_m49_country_regions.csv` | [UN M49](https://unstats.un.org/unsd/methodology/m49/) regions, including the nine-group classification used in Fig. 2D |
| `raw/controls/unhcr_refugee_stock_2022.csv` | [UNHCR](https://www.unhcr.org/refugee-statistics/) end-2022 refugee stock by country of asylum |
| `raw/controls/wdi_2019_2022_long.csv` | [World Development Indicators](https://databank.worldbank.org/source/world-development-indicators), 2019–2022, as retrieved for the paper |
| `raw/controls/wdi_employment_to_population_ratio_country_year.csv` | WDI `SL.EMP.TOTL.SP.ZS`, Gallup wave years |
| `raw/population/gadm1_ghs_pop_2025_population.csv` | [GHS-POP R2023A](https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php) epoch 2025 population summed within GADM 4.1 ADM1 units |
| `raw/geography/gadm1_subset_ids.csv` | GADM 4.1 ADM1 identifiers used by the pipeline (see below) |

### Obtain separately

| Path | Source and downloaded name | SHA-256 |
|---|---|---|
| `restricted/raw/Gallup_World_Poll_022026_ALL_WAVES.zip` | Gallup World Poll, all waves, archive 022026 (licensed; downloaded as `Gallup_World_Poll_022026_ALL WAVES_DAT_FILE.zip`) | `74588c6ff96a743c254430512d84d4dd9d42f3aa1d603e9c0c5316aaf992fc3c` |
| `restricted/raw/GWP_022026_Codebook_dta.txt` | Gallup codebook, archive 022026 (licensed) | `a1f3ed0ac16300c0797e390c7d6ccd14d82057b187014b2ae9ebe009c4d71f69` |
| `restricted/raw/GWP_022026_SysfileInfo_sav.xlsx` | Gallup variable and value-label metadata, archive 022026 (licensed) | `315380d874d7daf6462d0821f154f71413ed08677b1d38385b11c16aca137557` |
| `raw/meta/country_2026.csv` | [Meta SCI on HDX](https://data.humdata.org/dataset/social-connectedness-index), 2026 country-to-country (`country.csv`) | `0d9e13d555f281a2cb5a8afd3b2598a77af69aa90908509959c580bd0b74d5bf` |
| `raw/meta/all_region_to_country.zip` | Meta SCI on HDX, 2026 region-to-country (member `gadm1_to_country.csv`) | `44b99b667120f01b3d50bc0a3fa2a655c2ee3c1d1d9ac8488ba16e58e99c7893` |
| `raw/meta/country_country_aug2020.tsv` | August 2020 country SCI, file `Raw_data/SCI/country_country_aug2020.tsv` in the [replication data](https://data.mendeley.com/datasets/7wddm84w9r/1) of Bailey et al. (2021), "International trade and social connectedness" | `3a29f975baa30b25f9d6a80e4b6d1f878e958a3c822dbbb0576b0b9d248993cc` |
| `raw/geography/gadm_410_complete.gpkg` | [GADM 4.1](https://gadm.org/download_world.html), whole world as a single GeoPackage (`gadm_410-gpkg.zip`, containing `gadm_410.gpkg`) | `5a85ef31541c85e12eed7eabfebe59c88c05881367101696d86abef56bfb4ec2` |
| `raw/controls/undesa_pd_2024_ims.xlsx` | [UN DESA International Migrant Stock](https://www.un.org/development/desa/pd/content/international-migrant-stock), 2024 revision (`undesa_pd_2024_ims_stock_by_sex_destination_and_origin.xlsx`); the 2020 column is used | `0e10179d05186041a65cf5c6200943b2231701f9c1fd33cbbe21d23b8ea47316` |
| `raw/controls/Gravity_csv_V202211.zip` | [CEPII Gravity V202211](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8) | `7ffea0510a536587f2f9295bbafef4a312a36b366acc256eb6d9770fe35e44b1` |

**Authors' crosswalks.** Two hand-curated crosswalks from Gallup region codes to GADM
units are required by stages 01, 04, 07, and 08. They contain Gallup codebook region
labels, so they are not included; researchers with Gallup access can request them from
the authors.

| Path | Contents | SHA-256 |
|---|---|---|
| `raw/geography/v1/gallup_to_gadm1_crosswalk_v1.csv` | Gallup region code to GADM1 | `f4acf68f64d0253fc1709aeee2db8010767213e5053d413cfd54d4d45e2330a1` |
| `raw/geography/v2/prior_submarine_geography_crosswalk.csv` | Gallup region code to GADM1/GADM2 compositions | `cb4b82d48a3f2fe87184a5ed52af503dfa73120c4c0b15017ccef5ab638ef952` |

### Built by the pipeline

GADM's license does not allow redistributing its boundaries. Stage 00 builds
`raw/geography/gadm1_metadata.gpkg` from `gadm_410_complete.gpkg` by dissolving it to
the ADM1 units listed in `gadm1_subset_ids.csv` (about 10 minutes, first run only).
Stage 07 writes `raw/geography/v2/authoritative_composition_overrides.csv`.

## Conventions

- MAI is kept on its original 0–9 scale. GlobalSCI is the equal-weighted mean of log SCI
  to 177 foreign countries in the 2026 release, standardized over all 178 origins.
- Gallup weights are normalized to equal country totals within each wave; pooled models
  also give each wave equal total weight. Weights enter `lmer` as prior weights.
- Mixed models are fitted by maximum likelihood with Satterthwaite intervals. OLS
  sensitivities use country-clustered CR1 standard errors.
- All estimates are cross-sectional associations. The 2026 SCI release postdates every
  Gallup wave.
