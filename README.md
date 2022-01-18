# PUMP package

Please see the vignettes for examples of how to use this package.

Last updated: January 6, 2022.

Authors:
* Zarni Htet
* Kristen Hunter
* Luke Miratrix
* Kristin Porter

For randomized controlled trials (RCTs) with a single intervention being measured on multiple outcomes, researchers often apply a multiple testing procedure (such as Bonferroni or Benjamini-Hochberg) to adjust $p$-values.
Such an adjustment reduces the likelihood of spurious findings, but also changes the statistical power, sometimes substantially, which reduces the probability of detecting effects when they do exist.
However, this consideration is frequently ignored in typical power analyses, as existing tools do not easily accommodate the use of multiple testing procedures.
We introduce the PUMP R package as a tool for analysts to estimate statistical power, minimum detectable effect size, and sample size requirements for multi-level RCTs with multiple outcomes.
Multiple outcomes are accounted for in two ways.
First, power estimates from PUMP properly account for the adjustment in $p$-values from applying a multiple testing procedure.
Second, as researchers change their focus from one outcome to multiple outcomes, different definitions of statistical power emerge.
PUMP allows researchers to consider a variety of definitions of power, as some may be more appropriate for the goals of their study.
The package estimates power for frequentist multi-level mixed effects models, and supports a variety of commonly-used RCT designs and models and multiple testing procedures.
In addition to the main functionality of estimating power, minimum detectable effect size, and sample size requirements, the package allows the user to easily explore sensitivity of these quantities to changes in underlying assumptions.

For a full package description, including a detailed technical appendix, see https://arxiv.org/abs/2112.15273.
See also the shiny app https://mdrc.shinyapps.io/pump/.

We provide below one example of using PUMP to calculate a minimium detectable effect size (MDES).
The user specifies the RCT design and model (d_m), the multiple testing procedure (MTP, in this case Holm),
the target power (0.8), and the type of power desired (individual power for outcome 1).
The user also specifies a variety of design and model parameters, such as the number of outcomes, sample sizes at different levels, variation explained by covariates, etc.

```
m <- pump_mdes(
  d_m = "d3.2_m3fc2rc",         # choice of design and analysis strategy
  MTP = "HO",                   # multiple testing procedure
  target.power = 0.80,          # desired power
  power.definition = "D1indiv", # power type
  M = 5,                        # number of outcomes
  J = 3,                        # number of schools per block
  K = 21,                       # number districts
  nbar = 258,                   # average number of students per school
  Tbar = 0.50,                  # prop treated
  alpha = 0.05,                 # significance level
  numCovar.1 = 5,               # number of covariates at level 1
  numCovar.2 = 3,               # number of covariates at level 2
  R2.1 = 0.1, R2.2 = 0.7,       # explanatory power of covariates for each level
  ICC.2 = 0.05, ICC.3 = 0.4,    # intraclass correlation coefficients
  rho = 0.4 )                   # how correlated outcomes are
```