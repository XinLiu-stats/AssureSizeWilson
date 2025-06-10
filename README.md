# Sample size planning for proportions based on Wilson score confidence intervals with precision and assurance
This function estimates the required sample size to ensure that the lower bound of a confidence interval (CI) for either a single proportion or the difference between two independent proportions is no less than a specified threshold. 
It allows users to specify a desired assurance level and supports both the Wald, Wilson score, and Exact methods. Additionally, the function evaluates empirical assurance probabilities and empirical coverage probabilities through simulation studies to validate the results.

## Setup
Use the following command in R to install the package (which is the latest version):
```
library(devtools)
install_github("XinLiu-stats/AssureSizeWilson") 
```

## Sample szie estimation for a single proportion

### Usage
```
est_onesample(p, p0, alpha = 0.05, assurance = 0.8, rep = 10000)
```
`p`: Numeric. The true population proportion to be estimated. Must be in the range (0, 1). The value of p is usually determined based on pre-experimental or historical data.

`p0`: Numeric. The pre-specified threshold proportion. Must be in the range (0, 1). The pre-specified threshold proportion p0 is typically determined based on regulatory guidelines, industry standards, expert consensus, or historical data from similar studies.

`alpha`: Numeric. The significance level for the CI (e.g., 0.05 for 95\% confidence). Must be in the range (0, 1). Common choices are 0.05 or 0.01, depending on the desired level of confidence and regulatory requirements.

`assurance`: Numeric. The desired assurance probability (e.g., 0.8 for 80\% assurance). Must be in the range (0, 1). Typical values range from 0.80 to 0.95, depending on the level of assurance desired for achieving the study objective.

`rep`: Integer. The number of simulation replications to evaluate assurance and coverage probabilities. Must be a positive integer. The default is 10,000; larger values yield more precise estimates.

### Value

#### Input

`p`: The input `p`.	

`p0`: The input `p0`.

`alpha`: The input `alpha`.

`assurance`: The input `assurance`.	

`rep`: The input `rep`.

#### Output

`Nw`: Sample size estimated using the Wald method.

`Ns`: Sample size estimated using the Wilson score method.

`Ne`: Sample size estimated using the Exact method.

`EAPw`: Empirical assurance probability for the Wald method (percentage).

`EAPs`: Empirical assurance probability for the Wilson score method (percentage).

`EAPe`: Empirical assurance probability for the Exact method (percentage).

`ECPw`: Empirical coverage probability for the Wald method (percentage).

`ECPs`: Empirical coverage probability for the Wilson score method (percentage).

`ECPe`: Empirical coverage probability for the Exact method (percentage).

### Example
```
library(AssureSizeWilson) # load the AssureSizeWilson package

# estimate sample size for a single proportion
est_onesample(
  p = 0.937,
  p0 = 0.88,
  alpha = 0.05,
  assurance = 0.8,
  rep = 10000 )

# results
----------------------------------------------------------------------------------------
  p0     p    assurance   Nw    Ns   Ne    EAPw    EAPs    EAPe    CPw     CPs     CPe
 0.88  0.937    0.8      143   223   227   70.28   83.10   81.26  93.69   96.14   96.15
----------------------------------------------------------------------------------------

```

## Sample szie estimation for the difference between two independent proportions

### Usage
```
est_twosample(p1, p2, delta, r = 1, alpha = 0.05, assurance = 0.8, rep = 10000)
```
`p1`: Numeric. Numeric. True proportion in group 1 (e.g., experimental group). Must be in the range (0, 1). The value of p1 is usually determined based on pre-experimental or historical data.

`p2`: Numeric. True proportion in group 2 (e.g., control group). Must be in the range (0, 1). The value of p2 is usually determined based on pre-experimental or historical data.

`delta`: Numeric. The specified threshold for the lower limit of the confidence interval for (p1-p2). 

 To ensure clinical relevance and regulatory compliance, the selection of δ should follow the guidance below:
 
 ### General method
 Combine clinical judgment, historical data, and statistical reasoning, tailored to the primary endpoint. For Non-Inferiority Trials, the δ based on the standard treatment’s effect from historical placebo-controlled trials (M₁),  
  ensuring δ does not exceed M₁, and adjust to a more conservative margin (M₂) via clinical judgment.
 ### Regulatory Guidance
  When available, follow regulatory recommendations. Example: FDA proposed some noninferiority margins for some clinical endpoints (binary responses) such as cure rate for anti-infective drug products (e.g., topical antifungals or vaginal antifungals).
 ### No Prior Data
 In the absence of prior evidence, a standard effect size (i.e., effect size adjusted for standard deviation) between 0.25 and 0.5 is usually chosen as δ.

`r`: Numeric. Sample size allocation ratio (n1 / n2) between group 1 and group 2.

`alpha`: Numeric. The significance level for the CI (e.g., 0.05 for 95\% confidence). Must be in the range (0, 1). Common choices are 0.05 or 0.01, depending on the desired level of confidence and regulatory requirements.

`assurance`: Numeric. The desired assurance probability (e.g., 0.8 for 80\% assurance). Must be in the range (0, 1). Typical values range from 0.80 to 0.95, depending on the level of assurance desired for achieving the study objective.

`rep`: Integer. The number of simulation replications to evaluate assurance and coverage probabilities. Must be a positive integer. The default is 10,000; larger values yield more precise estimates.

### Value

#### Input

`p1`: The input `p1`.

`p2`: The input `p2`.

`delta`: The input `delta`.

`r`: The input `r`.

`alpha`: The input `alpha`.

`assurance`: The input `assurance`.

`rep`: The input `rep`.

#### Output

`n1_w`: Required sample sizes for group 1 using the Wald method.

`n2_w`: Required sample sizes for group 2 using the Wald method.

`n1_s`: Required sample sizes for group 1 using the Wilson score method.

`n2_s`: Required sample sizes for group 2 using the Wilson score method.

`n1_e`: Required sample sizes for group 1 using the Exact method.

`n2_e`: Required sample sizes for group 2 using the Exact method.

`EAPw`: Empirical assurance probability for the Wald method (percentage).

`EAPs`: Empirical assurance probability for the Wilson score method (percentage).

`EAPe`: Empirical assurance probability for the Exact method (percentage).

`ECPw`: Empirical coverage probability for the Wald method (percentage).

`ECPs`: Empirical coverage probability for the Wilson score method (percentage).

`ECPe`: Empirical coverage probability for the Exact method (percentage).

### Example
```
library(AssureSizeWilson) # load the AssureSizeWilson package

# estimate sample size for the difference between two independent proportions
est_twosample(
   p1 = 0.97,
   p2 = 0.98,
   delta = -0.06,
   r = 1,
   alpha = 0.05,
   assurance = 0.97,
   rep = 10000)

# results
------------------------------------------------------------------------------------------------------------------------
r   delta     p1     p2   assurance  n1_w  n2_w  n1_s  n2_s  n1_e  n2_e   EAPw    ECPw    EAPs   ECPs    EAPe   ECPe
1   -0.06    0.97   0.98   0.97      288   288   356   356   375   375    95.76   94.65   96.9   95.74   97.28  96.91
------------------------------------------------------------------------------------------------------------------------

