
# spatconf-functions

Functions for Semiparametric Approaches to Mitigate Spatial Confounding

Included are the functions to implement the simulation run in 'Semiparametric approaches for mitigating spatial confounding in large environmental epidemiology cohort studies'. 

| Function | Description |
|:----|:-------|
|simulation_script| Script used to run simulation |
|data_generating_functions| Functions used to create data |
|unadjusted_mod| Spatially-unadjusted model, no spatial splines used |
|gSEM_mod| Geoadditive Structural Equation method (Thaden and Kneib, 2018) |
|ks_mod| Method introduced by Keller and Szpiro (2020) |
|eps_mod| Exposure-Penalized Spline method (Bobb et al., 2022) |
|spatplus_mod| Spatial+ method (Dupont et al., 2022)  |
|dfspatplus_mod| Novel approach introduced in manuscript |

# References

Bobb, J. F., Cruz, M. F., Mooney, S. J., Drewnowski, A., Arterburn, D., and Cook, A. J.
(2022). Accounting for Spatial Confounding in Epidemiological Studies with Individual-
Level Exposures: An Exposure-Penalized Spline Approach. *Journal of the Royal Statistical
Society Series A: Statistics in Society*, 185(3):1271–1293.

Dupont, E., Wood, S. N., and Augustin, N. H. (2022). Spatial+: A novel approach to spatial
confounding. *Biometrics*, 78(4):1279–1290.

Keller, J. P. and Szpiro, A. A. (2020). Selecting a scale for spatial confounding adjustment.
*Journal of the Royal Statistical Society: Series A (Statistics in Society)*, 183(3):1121–1143.

Thaden, H. and Kneib, T. (2018). Structural Equation Models for Dealing With Spatial
Confounding. *The American Statistician*, 72(3):239–252.

