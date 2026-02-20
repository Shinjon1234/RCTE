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

## Monte Carlo Simulation

The simulation framework evaluates:

- ATE (Plug-in and Doubly Robust)
- ATT (Plug-in and Doubly Robust)
- ATO (Plug-in and Doubly Robust)
- RCTE (Plug-in and Doubly Robust)

Across 16 combinations of:
- Treatment Effect specification
- Propensity score specification

File:
simulations/mc_full_simulation.R

## Empirical Application: LaLonde Data

We analyze the NSW job training program using:

- RCTE (Plug-in and Doubly Robust)
- ATO
- ATE (AIPW)
- ATT (Augmented)

File:
application/lalonde_analysis.R

## Empirical Application: 401(k) Eligibility

We analyze the effect of 401(k) eligibility on net financial assets using:

- RCTE (Plug-in and Doubly Robust)
- ATO
- ATE (AIPW)
- ATT (AIPW)

File:
application/401k_analysis.R




