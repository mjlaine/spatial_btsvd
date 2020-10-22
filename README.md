# Matlab files for Bayesian truncated SVD analysis with spatial effects

Matlab code for analyses in the article:
V. Junttila, M. Laine: Bayesian principal component regression model with spatial effects for forest inventory under small field sample size, *Remote Sensing of Environment*, **192**,
pages 45-57, 2017. [doi:10.1016/j.rse.2017.01.035](https://doi.org/10.1016/j.rse.2017.01.035)

The code estimates a linear random effect model y = Xβ + η + ε, where we have correlated zero-mean random effect η ~ N(0,C) and uncorrelated noise ε ~ N(0,τ²I). Prior for the regression parameter β is β ∼ N(0, αᵢ¯¹), i = 0,1,...,p, with hyper prior for precision parameters α as αᵢ ∼ χ²(νᵢ, aᵢ), i = 0,1,...,p, where χ²(νᵢ, aᵢ) is the scaled chi squared distribution. Spatial covariance matrix C is defined by correlation decay parameter φ and spatial variance σ², thus making the total parameter vector for this model as θ = [β,α,τ,φ,σ]. For the regression parameters, principle component regression by truncated SVD is performed. See the article for more details.

## Matlab code

The main functions provided here are

```
[results,chain,sschain] = spatial_mcmcrun(x,y,xcoord,ycoord,opts)
```

for MCMC analysis of the model parameters,

```
out=spatial_prediction(results,chain,xnew,xcoordnew,ycoordnew,nsimu)
```

for predicting the output at new spatial locations,

```
spatial_mcmcplot(chain,results,type,varargin)
```

for plotting various MCMC chain plots, and

```
[w,a,b]=BtSVD2(x,y,v)
```

for Bayesian truncated SVD regression.


## Test case with synthetic data

You can test the code with the following two scripts:

`makesyntdata.m`: generate synthetic spatial data set.

`analyzesynt.m`: calculate predictions by different methods.

The test code uses some auxiliary Matlab functions. Some statistical distributions from the MCMCSTAT toolbox, which is available from GitHub at
<https://www.github.com/mjlaine/mcmcstat>, and
some linear regression routines from the RSTOOLS toolbox, available at
<https://www.github.com/mjlaine/rstools> .

