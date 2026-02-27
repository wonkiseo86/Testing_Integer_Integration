# R code for "Testing for integer integration in functional time series"

---
# Setup Instructions for Reproducibility

To ensure full reproducibility of the results presented in the paper, please configure your environment as follows:
### Software Requirements
The analysis was performed using **R version 4.5.1**. The following R packages and their specific versions are required to replicate the results:

* **Packages:** `fda` (6.3.0), `tseries` (0.10-58), `sandwich` (3.1-1), `sde` (2.0.18), `variables` (1.1-2), `basefun` (1.2-3), `polynom` (1.4-1), `fracdiff` (1.5-3), `arfima` (1.8-1), `fdaACF` (1.0.0), `ftsa` (6.6), and `LaplacesDemon` (16.1.6).

### Data Preparation
All data files must be placed in a directory named `data` at the project root.
* **Step 1:** Create a folder named `data` in the same directory as this README.
* **Step 2:** Locate the provided `data.zip` file.
* **Step 3:** Unzip the contents of `data.zip` directly into the `data` folder.
  
> **Note:** After unzipping, the file paths should look like this:
> - `data/Canadian_daily_yields.csv`
> - `data/NUT2/1f.txt`, `data/NUT2/1m.txt`, etc.

## 📂 Files and Scripts

### Directories & Archives
* **auxiliary**: Sub-folder containing some generic functions
* **data.zip**: Folder containing Canadian yield curve data and French mortality data

### Empirical Study Scripts
* **empirical_canadian_yield.r**: Study on the integration property of Canadian yield curves
* **empirical_French_mortality_revision.r**: Study on the integration property of French age-specific mortality rates
* **empirical_French_mortality_Supplement_revision.r**: Supplementary results on French age-specific mortality rates
* **figures.r**: Create Figures 1 and 2
  
### Simulation Study Scripts
* **Int_integ_sim_accuracy_fraction_revision.r**: Accuracy of the testing procedure when the DGP is fractionally-integrated
* **Int_integ_sim_accuracy_integer_revision.r**: Accuracy of the testing procedure when the DGP is integer-integrated
* **Int_integ_sim_size_power_properties_revision.r**: Size and power of the tests
* **Int_integ_sim_size_power_properties_bandwidths.r**: Size and power of the tests depending on the bandwidth parameter choice
* **Int_integ_sim_size_power_properties_demeaned_revision.r**: Size and power of the tests computed with the demeaned time series
* **Int_integ_sim_size_power_properties_local_power_d_revision**: Power against local alternatives

### Data Files
* **CV1.Rdata**: Critical values for various quantiles which are used for the Canadian yield curve example (to compute approximate p-values)
* **Critical_values.r**: R script used to compute Critical values for various quantiles which are used for the Canadian yield curve example. 

---

## 📖 Detailed Instruction for Reproduction of Results

### 📊 Empirical Results

#### 1. Table 4 of the Main Article
* **Data Required:** `Canadian_daily_yields.csv`
* **Script:** `empirical_canadian_yield.r`
* **Description:** Results are directly obtained from the script.

#### 2. Table 5 of the Main Article
* **Data Required:** Datasets in the `NUT2` folder
* **Scripts:** `empirical_French_mortality_revision.r` and `empirical_French_mortality_supplement_revision.r`
* **Description:**
    * (i) Results for the subregions except the 9th region (Alsace) are obtained using `empirical_French_mortality_revision.r`. Set the `transformation` parameter as `transformation=1`.
    * (ii) Results for the 9th region (Alsace) are obtained using Section 1 of `empirical_French_mortality_supplement_revision.r`. Set the `transformation` parameter as `transformation=1`.

#### 3. Table 1 of the Supplementary Material
* **Data Required:** Datasets in the `NUT2` folder
* **Scripts:** `empirical_French_mortality_revision.r` and `empirical_French_mortality_supplement_revision.r`
* **Description:**
    * (i) Results for the subregions except the 9th region (Alsace) are obtained using `empirical_French_mortality_revision.r`. Set the `transformation` parameter as `transformation=3`.
    * (ii) Results for the 9th region (Alsace) are obtained using Section 1 of `empirical_French_mortality_supplement_revision.r`. Set the `transformation` parameter as `transformation=3`.

#### 4. Table 2 of the Supplementary Material
* **Data Required:** Datasets in the `NUT2` folder
* **Scripts:** `empirical_French_mortality_revision.r` and `empirical_French_mortality_supplement_revision.r`
* **Description:**
    * (i) Results for the subregions except the 9th region (Alsace) are obtained using `empirical_French_mortality_revision.r`. Set the `transformation` parameter as `transformation=4`.
    * (ii) Results for the 9th region (Alsace) are obtained using Section 1 of `empirical_French_mortality_supplement_revision.r`. Set the `transformation` parameter as `transformation=4`.
 
#### 5. Figures 1 and 2
* **Data Required:** `Canadian_daily_yields.csv` and Datasets in the `NUT2` folder
* **Scripts:** `figures.r`
* **Description:**
    * Figures 1 and 2 are generated from the code.
      
---

### 🧪 Simulation Results

#### 1. Table 1 of the Main Article
* **Script:** `Int_integ_sim_size_power_properties_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and for each test (V0 test and V1 test).
    * (ii) Results are obtained using `Int_integ_sim_size_power_properties_revision.r`. Set the `bdd` parameter (corresponding to the value b) and the `nonstat` parameter (determining V0 or V1 test) accordingly. For example, `bdd=0.6` and `nonstat=0` (results for V0 test when b=0.6); `bdd=0.15` and `nonstat=1` (results for V1 test when b=0.15).

#### 2. Table 2 of the Main Article
* **Scripts:** `Int_integ_sim_size_accuracy_integer_revision.r` and `Int_integ_sim_size_fraction_integer_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6).
    * (ii) Results for Table 2(a) are obtained using `Int_integ_sim_size_accuracy_integer_revision.r`. Set the `bdd` parameter (corresponding to b) as `bdd=0.15` for upper rows and `bdd=0.6` for lower rows.
    * (iii) Results for Table 2(b) are obtained in the same way but using `Int_integ_sim_size_fraction_integer_revision.r`.

#### 3. Table 3 of the Main Article
* **Script:** `Int_integ_sim_size_power_properties_demeaned_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6).
    * (ii) Results for Table 3 are obtained using `Int_integ_sim_size_power_properties_demeaned_revision.r`. Set `bdd=0.15` for upper rows and `bdd=0.6` for lower rows.

#### 4. Table 3 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_bandwidths.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and the bandwidth q.
    * (ii) First, set `uband=0` (determining bandwidth q). Then, set `bdd=0.15` for upper rows and `bdd=0.6` for lower rows.

#### 5. Table 4 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_bandwidths.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of b (= 0.15, 0.6) and the bandwidth q.
    * (ii) First, set `uband=1/4`. Then, set `bdd=0.15` for upper rows and `bdd=0.6` for lower rows.

#### 6. Table 5 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_bandwidths.r`
* **Description:**
    * (i) **Note:** Results across various values of T and d_1 are obtained together for each value of the bandwidth q.
    * (ii) Set `bdd=0.75`. Then, set `uband=0`, `uband=1/5` and `uband=1/4` to obtain results for cases where q ~ log T, q ~ T^{1/5}, and q ~ T^{1/4}, respectively.

#### 7. Table 6 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_local_power_d_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).
    * (ii) Set `uband=0`. Then, set `bdd=0.15` (resp. `bdd=0.6`) for upper (resp. lower) rows.

#### 8. Table 7 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_local_power_d_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).
    * (ii) Set `uband=1/5`. Then, set `bdd=0.15` (resp. `bdd=0.6`) for upper (resp. lower) rows.

#### 9. Table 8 of the Supplementary Material
* **Script:** `Int_integ_sim_size_power_properties_local_power_d_revision.r`
* **Description:**
    * (i) **Note:** Results across various values of T and c are obtained together for each value of b (= 0.15, 0.6).
    * (ii) Set `uband=1/4`. Then, set `bdd=0.15` (resp. `bdd=0.6`) for upper (resp. lower) rows.

#### 10. Figure 1 of the Supplementary Material
* Produced from the results for Tables 6-8.
