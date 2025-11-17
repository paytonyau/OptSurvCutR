# OptSurvCutR: Optimal Survival Cut-Point Discovery for Time-to-Event Analysis

Provides a robust 3-step workflow for optimal cut-point analysis in
time-to-event (survival) data. Functions determine the optimal number of
cut-points ('find_cutpoint_number()'), find their precise locations
('find_cutpoint()') using systematic or genetic algorithms (via
'rgenoud'), and validate stability via bootstrapping
('validate_cutpoint()'). Analyses can be adjusted for covariates using
standard 'survival' package models. Ideal for biomarker analysis and
patient stratification, where non-linear relationships may exist.
Methods described in Yau, P.T.O. (2025)
[doi:10.1101/2025.10.08.681246](https://doi.org/10.1101/2025.10.08.681246)
.

## See also

Useful links:

- <https://github.com/paytonyau/OptSurvCutR>

- <https://paytonyau.github.io/OptSurvCutR/>

- Report bugs at <https://github.com/paytonyau/OptSurvCutR/issues>

## Author

**Maintainer**: Payton Yau <tungon@gmail.com>
([ORCID](https://orcid.org/0000-0002-3283-0370)) \[copyright holder\]
