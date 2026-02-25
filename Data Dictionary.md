# Data Dictionary

This repository contains two primary datasets used for the replication study. Below are the detailed definitions for each variable and the nature of the observations.


## 1. Yield Curve Dataset
**Definition of Observation:** Each row represents a daily snapshot of the zero-coupon yield curve based on market closing prices for a specific business day.

| Variable Name | Data Type | Description |
| :--- | :--- | :--- |
| **Date** | Date | The specific trading date of the observation (MM/DD/YYYY). |
| **ZC025YR** | Numeric | 0.25-year (3-month) zero-coupon bond yield (decimal form). |
| **ZC050YR** | Numeric | 0.5-year (6-month) zero-coupon bond yield (decimal form). |
| **ZC100YR** | Numeric | 1-year zero-coupon bond yield (decimal form). |
| **ZCXXXYR** | Numeric | Zero-coupon yield for the maturity of XXX years (e.g., ZC650YR = 65 years). |

---

## 2. Mortality Dataset (French Regions)
**Definition of Observation:** Each row represents annual mortality indicators for a specific age and year within a demographic group defined by region and gender.

### File Naming Convention
- **Numeric Prefix (1–22):** Indicates the **Region Code** (see mapping below).
- **Suffix (m/f):** Denotes gender (**m** for Male, **f** for Female).
- *Example: `1f.csv` contains mortality data for females in Île de France.*

### Region Code Mapping
| Code | Region Name | Code | Region Name |
| :--- | :--- | :--- | :--- |
| 1 | Île de France | 12 | Pays de la Loire |
| 2 | Centre-Val de Loire | 13 | Bretagne |
| 3 | Bourgogne | 14 | Aquitaine |
| 4 | Franche-Comté | 15 | Limousin |
| 5 | Basse-Normandie | 16 | Poitou-Charentes |
| 6 | Haute-Normandie | 17 | Languedoc-Roussillon |
| 7 | Nord-Pas-de-Calais | 18 | Midi-Pyrénées |
| 8 | Picardie | 19 | Auvergne |
| 9 | Alsace | 20 | Rhône-Alpes |
| 10 | Champagne-Ardenne | 21 | Provence-Alpes-Côte d’Azur |
| 11 | Lorraine | 22 | Corse |

### Variable Definitions
| Variable | Description |
| :--- | :--- |
| **Year** | The calendar year of the observation. |
| **Age** | The exact age $x$ at the beginning of the 1-year age interval. |
| **mx** | Central death rate between ages $x$ and $x+1$. |
| **qx** | Probability of death between ages $x$ and $x+1$. |
| **ax** | Average length of survival for persons dying in the interval. |
| **lx** | Number of survivors at exact age $x$ (standardized to a radix of 100,000). |
| **dx** | Number of deaths between ages $x$ and $x+1$. |
| **Lx** | Number of person-years lived between ages $x$ and $x+1$. |
| **Tx** | Total number of person-years lived after exact age $x$. |
| **ex** | Life expectancy at exact age $x$ (in years). |
