# eiu-wellbeing-replication

Replication code for: Zusammenhang zwischen den fünf Dimensionen des EIU Democracy Index (EIU-DI) und der internationalen Lebenszufriedenheit (Grewing 2026, University of Oldenburg).

## Files
- `FINAL_R_SKRIPT_Grewing_Felix.R` — main analysis script (OLS main models, V-Dem robustness check, diagnostics)
- `sessionInfo.txt` — R 4.5.1 package versions

## Data

Data files are not included due to licensing. Place them in the working directory before running the script.

Main data:
- EIU Democracy Index 2024 — proprietary (eiu.com)
- Life Satisfaction, Social Support — [World Happiness Report 2025] (https://worldhappiness.report/data-sharing/)
- V-Dem indices v15 — [V-Dem Dataset Archive → Core v15] (https://www.v-dem.net/data/dataset-archive/)

Control variables:
Control variables:
- ln(GDP per capita, PPP) — [World Bank] (https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.KD)
- Unemployment rate — [World Bank] (https://data.worldbank.org/indicator/SL.UEM.TOTL.ZS)
- Age dependency ratio — [World Bank] (https://data.worldbank.org/indicator/SP.POP.DPND.OL)
- Social protection expenditure (% GDP) — [QoG Standard Dataset Jan26] (https://www.gu.se/en/quality-government/qog-data/data-downloads/standard-dataset) (originally ILO)
- WGI Aggregate (6 indicators) — [QoG Standard Dataset Jan26] (https://www.gu.se/en/quality-government/qog-data/data-downloads/standard-dataset) (originally World Bank)
- Social Support — [World Happiness Report 2025] (https://worldhappiness.report/data-sharing/)

V-Dem v16 was released March 2026; this analysis uses v15.

## How to run

Open the script in R (≥ 4.5.1), place data files in the working directory, run from top to bottom. Results go to `out_analyse_<timestamp>/out/`. Reproducible via `set.seed(123)`.
