# Risk Calibrated Treatment Effect (RCTE):
This repository contains the full theoretical derivations, simulation studies,
and empirical applications for the weighted average treatment effect (WATE)
with weight:

    w(X) = 1 − μ0(X)

The estimand studied is:

    τ_w = E[w(X)(Y(1) − Y(0))]

We construct a plug-in estimator and a doubly robust estimator, derive its efficient influence
function (EIF), and study its asymptotic properties under semiparametric efficiency theory. 
We have also studied the small sample properties of the estimators and establish their small sample efficiency 
using concentration bounbds from ststistical learning theory.

---

## Methodological Contributions

• Formal definition of RCTE-type WATE
• Doubly robust estimator construction
• Derivation of semiparametric efficiency
• Asymptotic normality proof
• Concentration bounds derivation

---

## Repository Structure

R/
  Core reusable functions:
    - Data generating processes
    - Nuisance model estimation
    - Weight construction
    - DR estimator
   

simulations/
  Monte Carlo experiments and bootstrap studies

application/
  Real data empirical analysis

manuscript/
  RMarkdown source for the paper and supplement

output/
  Generated tables and figures

---

rcte_dr_wate/
│
├── simulations/
│   ├── mc_design_reference.R
│   ├── mc_run_all_scenarios.R
│   └── mc_helpers.R




