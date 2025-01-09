# Sample size planning for proportions based on Wilson score confidence intervals with precision and assurance
This function estimates the required sample size to ensure that the lower bound of a confidence interval (CI) for either a single proportion or the difference between two independent proportions is no less than a specified threshold. 
It allows users to specify a desired assurance level and supports both the Wald and Wilson score methods. Additionally, the function evaluates empirical assurance probabilities and CI coverage probabilities through simulation studies to validate the results.

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
`p`: Numeric. The true population proportion to be estimated. Must be in the range (0, 1).

`p0`: Numeric. The pre-specified threshold proportion. Must be in the range (0, 1).

`alpha`: Numeric. The significance level for the CI (e.g., 0.05 for 95\% confidence). Must be in the range (0, 1).

`assurance`: Numeric. The desired assurance probability (e.g., 0.8 for 80\% assurance). Must be in the range (0, 1).

`rep`: Integer. The number of simulation replications to evaluate assurance and coverage probabilities. Must be a positive integer.

### Value
`p0`: The input `p0`.

`p`: The input `p`.	

`assurance`: The input `assurance`.	

`Nw`: Sample size estimated using the Wald method.

`Ns`: Sample size estimated using the Wilson score method.

`EAPw`: Empirical assurance probability for the Wald method (percentage).

`EAPs`: Empirical assurance probability for the Wilson score method (percentage).

`CPw`: Confidence interval coverage probability for the Wald method (percentage).

`CPs`: Confidence interval coverage probability for the Wilson score method (percentage).

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
------------------------------------------------------------------
  p0     p    assurance   Nw    Ns    EAPw    EAPs   CPw     CPs
 0.88  0.937    0.8      143   223   71.00   83.61  94.47   96.23
------------------------------------------------------------------

```

## Sample szie estimation for the difference between two independent proportions

### Usage
```
est_twosample(p1, p2, delta, r = 1, alpha = 0.05, assurance = 0.8, rep = 10000)
```
`p1`: Numeric. Numeric. True proportion in group 1 (e.g., experimental group). Must be in the range (0, 1).

`p2`: Numeric. True proportion in group 2 (e.g., control group). Must be in the range (0, 1).

`delta`: Numeric. The specified threshold for the lower limit of the confidence interval for (p1-p2).

`r`: Numeric. Sample size allocation ratio (n1 / n2) between group 1 and group 2.

`alpha`: Numeric. The significance level for the CI (e.g., 0.05 for 95\% confidence). Must be in the range (0, 1).

`assurance`: Numeric. The desired assurance probability (e.g., 0.8 for 80\% assurance). Must be in the range (0, 1).

`rep`: Integer. The number of simulation replications to evaluate assurance and coverage probabilities. Must be a positive integer.

### Value
`r`: Numeric. The input `r`.

`delta`: Numeric. The input `delta`.

`p1`: Numeric. The input `p1`.

`p2`: Numeric. The input `p2`.

`assurance`: Numeric. The input `assurance`.

`n1_w`: Required sample sizes for group 1 using the Wald method.

`n2_w`: Required sample sizes for group 2 using the Wald method.

`n1_s`: Required sample sizes for group 1 using the Wilson score method.

`n2_s`: Required sample sizes for group 2 using the Wilson score method.

`EAPw`: Empirical assurance probability for the Wald method (percentage).

`EAPs`: Empirical assurance probability for the Wilson score method (percentage).

`CPw`: Confidence interval coverage probability for the Wald method (percentage).

`CPs`: Confidence interval coverage probability for the Wilson score method (percentage).

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
--------------------------------------------------------------------------------------------------
r   delta     p1     p2   assurance  n1_w  n2_w  n1_s  n2_s   EAPw    ECPw   EAPs   ECPs
1   -0.06    0.97   0.98   0.97      288   288   356   356   95.76   94.65   96.9   95.74
--------------------------------------------------------------------------------------------------

