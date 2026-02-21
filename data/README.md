# Data Description

This folder documents the datasets used in the empirical applications of the paper.

No raw data is stored locally in this repository. All datasets are publicly available and are programmatically downloaded within the scripts to ensure reproducibility.

---

## 1. LaLonde Dataset

Source:
https://users.nber.org/~rdehejia/data/

Files used:
- nswre74_treated.txt
- cps_controls.txt

Description:
The LaLonde dataset is a benchmark dataset in causal inference evaluating the impact of a job training program (NSW) on earnings.

Variables:
- treat: Treatment indicator (1 = treated)
- age: Age
- educ: Years of education
- black: Indicator for Black ethnicity
- hisp: Indicator for Hispanic ethnicity
- married: Marital status
- nodegr: No high school degree
- re74: Earnings in 1974
- re75: Earnings in 1975
- re78: Earnings in 1978 (outcome)

Usage in this repository:
- Outcome: re78
- Treatment: treat
- Covariates: age, educ, black, hisp, married, nodegr, re74, re75

---

## 2. 401(k) Dataset

Source:
DoubleML R package

Accessed via:
fetch_401k()

Description:
This dataset studies the causal effect of 401(k) eligibility on household net financial assets.

Variables:
- e401: Eligibility for 401(k) (treatment)
- net_tfa: Net financial total assets (outcome)
- age: Age
- inc: Income
- educ: Education
- fsize: Family size
- marr: Marital status
- twoearn: Dual-income household indicator
- db: Defined benefit pension
- pira: IRA participation
- hown: Home ownership

Usage in this repository:
- Outcome: net_tfa
- Treatment: e401
- Covariates: age, inc, educ, fsize, marr, twoearn, db, pira, hown

---

## 3. Simulated Data

Simulated datasets are generated within the repository using:

scripts:
- simulations/mc_full_simulation.R

Data Generating Process (DGP):
- Covariates: Multivariate normal
- Treatment: Logistic / nonlinear propensity score
- Outcomes: Structured via μ₀(X) and τ(X)

No external data is required for simulations.

---

## Reproducibility Notes

All datasets are:
- Publicly available
- Automatically downloaded or generated
- Fully reproducible via provided scripts

No manual data download is required.
