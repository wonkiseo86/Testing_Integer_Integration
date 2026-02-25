R code for "Testing for integer integration in functional time series"
<br><br>
################ Files and Scripts ################

auxiliary: sub-folder containing some generic functions

data.zip: folder containing Canadian yield curve data and French mortality data

empirical_canadian_yield.r: study on the integration property of Canadian yield curves

empirical_French_mortality_revision.r: study on the integration property of French age-specific mortality rates

empirical_French_mortality_Supplement_revision.r: supplementary results on French age-specific mortality rates

Int_integ_sim_accuracy_fraction_revision.r: accuracy of the testing proceudre when the DGP is fractionally-integrated

Int_integ_sim_accuracy_integer_revision.r: accuracy of the testing proceudre when the DGP is integer-integrated

Int_integ_sim_size_power_properties_revision.r: size and power of the tests

Int_integ_sim_size_power_properties_bandwidths.r: size and power of the tests depending on the bandwidth parameter choice

Int_integ_sim_size_power_properties_demeaned_revision.r: size and power of the tests computed with the demeaned time series

Int_integ_sim_size_power_properties_local_power_d_revision: power against local alternatives

CV1.Rdata: critical values for various quantiles which are used for the Canadian yield curve example (to compute approximate the p-values). 

 
 
<br> <br><br>
######### Detailed Instruction for Reproduction of the Empirical/Simulation Results #########
<br><br>
Empirical results: <br>
1. Table 4 of the Main Article<br>
Data Required: "Canadian_daily_yields.csv".<br> 
Script: "empirical _canadian_yield.r".<br> 
Description: Results are directly obtained from the script.

2. Table 5 of the Main Article<br>
Data Required: Datasets in the "NUT2" folder. <br>
Scripts: "empirical _French_mortality_revision.r" and "empirical _French_mortality_supplement_revision.r"<br>
Description:<br>
(i) Results for the subregions except the 9th region (Alsace) are obtained using "empirical _French_mortality_revision.r". Set the "transformation" parameter  as "transformation=1". <br>
(ii) Results for the 9th region (Alsace) are obtained using Section 1 of "empirical _French_mortality_supplement_revision.r". Set the "transformation" parameter  as "transformation=1". <br>

3. Table 1 of the Supplementary Material<br>
Data Required: Datasets in the "NUT2" folder. <br>
Scripts: "empirical _French_mortality_revision.r" and "empirical _French_mortality_supplement_revision.r"<br>
Description:<br>
(i) Results for the subregions except the 9th region (Alsace) are obtained using "empirical _French_mortality_revision.r". Set the "transformation" parameter  as "transformation=3".<br>
(ii)  Results for the 9th region (Alsace) are obtained using Section 1 of "empirical _French_mortality_supplement_revision.r". Set the "transformation" parameter  as "transformation=3". <br>

4. Table 2 of the Supplementary Material<br>
Data Required: Datasets in the "NUT2" folder. <br>
Scripts: "empirical _French_mortality_revision.r" and "empirical _French_mortality_supplement_revision.r"<br>
Description:<br>
(i) Results for the subregions except the 9th region (Alsace) are obtained using "empirical _French_mortality_revision.r". Set the "transformation" parameter  as "transformation=4".<br>
(ii)  Results for the 9th region (Alsace) are obtained using Section 1 of "empirical _French_mortality_supplement_revision.r". Set the "transformation" parameter  as "transformation=4".<br>
<br><br>
Simulation results:<br>
1. Table 1 of the Main Article<br>
Script: "Int_integ_sim_size_power_properties_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and for each test (V0 test and V1 test).<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_revision.r". Set the  "bdd" parameter (corresponding to the value b) and the "nonstat" parameter (determining V0 or V1 test) accordingly. For example, "bdd=0.6" and "nonstat=0" (results for V0 test when b=0.6); "bdd=0.15" and "nonstat=1" (results for V1 test when b=0.15) <br>

2. Table 2 of the Main Article<br>
Scripts: "Int_integ_sim_size_accuracy_integer_revision.r" and "Int_integ_sim_size_fraction_integer_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6).<br>
(ii) Results for Table 2(a) are obtained using "Int_integ_sim_size_accuracy_integer_revision.r". Set the "bdd" parameter (corresponding to the value b) as "bdd=0.15" to obtain the upper rows of Table 2(a) and "bdd=0.6" to obtain the lower rows of Table 2(a). <br>
(iii) Results for Table 2(b) are obtained in the same way but using "Int_integ_sim_size_fraction_integer_revision.r"<br>

3. Table 3 of the Main Article<br>
Script: "Int_integ_sim_size_power_properties_demeaned_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6).<br>
(ii) Results for Table 3 are obtained using "Int_integ_sim_size_power_properties_demeaned_revision.r". Set the "bdd" parameter (corresponding to the value b) as "bdd=0.15" to obtain the upper rows and "bdd=0.6" to obtain the lower rows.  <br>

3. Table 3 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_bandwidths.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and the bandwidth q.<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_bandwidths.r". First, set the "uband" parameter (determining bandwidth q) as "uband=0". Then, set the "bdd" parameter   (corresponding to the value b) as "bdd=0.15" to obtain the upper rows and "bdd=0.6" to obtain the lower rows.  <br>

4. Table 4 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_bandwidths.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and the bandwidth q.<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_bandwidths.r". First, set the "uband" parameter (determining bandwidth q) as "uband=1/4". Then, set the "bdd" parameter  (corresponding to the value b) as "bdd=0.15" to obtain the upper rows and "bdd=0.6" to obtain the lower rows.  <br>

5. Table 5 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_bandwidths.r"<br>
Description:<br>
(i) Note: Results across various values of T and d_1 are obtained together for each value of the bandwidth q.<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_bandwidths.r". First set the "bdd" parameter (corresponding to the value b) as "bdd=0.75". Then, set the "uband" parameter (determining bandwidth q) as "uband=0", "uband=1/5" and "uband=1/4" to obtain the results for the cases where q ~ log T, q ~ T^{1/5}, and q ~ T^{1/4}, respectively.<br>

6. Table 6 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_local_power_d_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_local_power_d_revision.r". First set the "uband" parameter (determining bandwidth q) as "uband=0". Then, set the "bdd" parameter (corresponding to the value b) as "bdd=0.15" (resp. "bdd=0.6") to obtain the results for the upper (resp. lower) rows.<br>

7. Table 7 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_local_power_d_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_local_power_d_revision.r". First set the "uband" parameter (determining bandwidth q) as "uband=1/5". Then, set the "bdd" parameter (corresponding to the value b) as "bdd=0.15" (resp. "bdd=0.6") to obtain the results for the upper (resp. lower) rows.<br>

8. Table 8 of the Supplementary Material<br>
Script: "Int_integ_sim_size_power_properties_local_power_d_revision.r"<br>
Description:<br>
(i) Note: Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).<br>
(ii) Results are obtained using "Int_integ_sim_size_power_properties_local_power_d_revision.r". First set the "uband" parameter (determining bandwidth q) as "uband=1/4". Then, set the "bdd" parameter (corresponding to the value b) as "bdd=0.15" (resp. "bdd=0.6") to obtain the results for the upper (resp. lower) rows.<br>

9. Figure 1 of the Supplementary Material is produced from the results for Tables 6-8.  <br>





# Data Dictionary

This repository contains two primary datasets used for the replication study. Below are the detailed definitions for each variable and the nature of the observations.

---

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
