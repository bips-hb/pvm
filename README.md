## PharmacoVigilance Methods (PVM)

`PVM` is an `R` package containing a wide variety of methods used in the field of pharmacovigilance for discovering 'interesting' drug-adverse event pairs from spontaneous reporting data. 

### Available Methods

The methods currently implemented are: 

1. the reporting odds ratio (ROR), see `R/ROR.R`;

2. Yule's Q, see `R/YulesQ.R`;

3. the proportional relative risk (PRR), see `R/PRR.R`;

4. the relative report rate (RRR), see `R/RRR.R`;

5. the reporting Fisher's exact test (RFET) and the mid-p-value test (midRFET), see `R/fisherExactTest.R`;

6. the chi-squared test (with and without Yates' correction for continuity), see `R/chi2Test.R`; 

7. the binomial likelihood ratio test, see `R/logLikelihoodRatioBinomial.R`; 

8. the test of the Poisson mean, see `R/PoissonTest.R`;

9. the Bayesian confidence propagation neural network (BCPNN), see `R/BCPNN.R`;

10. the Gamma Poisson shrinker (GPS), see `R/GPS.R`, and 

11. the LASSO, see `R/LASSO.R`

### Installation 
To install, simply type in R

```R
devtools::install_github("bips-hb/pvm")
```

### References

Please cite 

__Adverse Event Discovery for Spontaneous Reporting Systems: A Systematic Comparison__\
*L.J. Dijkstra, M. Garling, R. Foraita & I. Pigeot*\
To be submitted (2018)

### Contact

Louis Dijkstra\
Leibniz Institute for Prevention Research & Epidemiology  
E-mail: dijkstra (at) leibniz-bips.de