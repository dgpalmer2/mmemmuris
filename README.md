# mmemmuris: Marginal and Mixed-Effects Models for Multivariate and Univariate Results Including Sphericity

## Pronunciation

The pronunciation of `mmemmuris` is "memories".

## Background

With the advent of linear mixed-effects models (LMM), multivariate analysis of variance (MANOVA) and repeated measures ANOVA (RM ANOVA) are falling out of popularity.  MANOVA and RM ANOVA have a few main drawbacks that deter statisticians from using it:

1) Listwise deletion for MANOVA/RM ANOVA data in the wide format causes an entire experimental unit with at least one missing value to be dropped from the analysis.  Marginal models and LMM process the data in the long format, so instead of dropping an entire experimental unit with missing data, only the time point(s) with missing data is dropped.

2) MANOVA only uses models with unstructured covariance structures (multivariate tests) and RM ANOVA only uses models with compound symmetry covariance structures (univariate tests).  Marginal models and LMM can fit models with these covariance structures and many others.

3) Multiple comparisons cannot be done on the between or within-subject factors unless the data is then reshaped from the wide format to the long format.  More often than not, researchers want to make inference on these factors.

4) MANOVA/RM ANOVA cannot treat the time variable as continuous, add polynomial terms for the time variable, or drop the interaction from the model, but marginal models and LMM can.

## Goals

The goals of this package are to:

1) Get researchers to start becoming familiar with marginal models and LMM through their relationships with MANOVA and RM ANOVA.

2) Get researchers to get comfortable enough with marginal models and LMM to use them as their new tools of choice instead of MANOVA and RM ANOVA.

## Prequisites

I'm assuming you've analyzed data using RM ANOVA and/or MANOVA before in R (`stats` and/or `car` packages) or some other software.

## Usage

To use this package:

1) Install either the `remotes` package or `devtools` package from CRAN.

        install.packages("remotes")
        install.packages("devtools")


2) Load either the `remotes` package or `devtools` package.

        library(remotes)
        library(devtools)

3) 

    a. Use the `install_github` function to install my package along with its dependencies (`nlme`, `MASS`, and `Matrix` packages). 
    
            install_github("dgpalmer2/mmemmuris")

    b. Alternatively, install the `knitr` and `rmarkdown` packages along with the above packages from CRAN to get access to the vignette guide on how to use my package.


            install_github("dgpalmer2/mmemmuris", dependencies = TRUE, build_vignettes = TRUE)
            browseVignettes("mmemmuris")
            
4) Load my package.

        library(mmemmuris)
        
## Documentation

A brief description of the functions used in my package along with their arguments and mathematical theory can be accessed with a question mark preceding the function name.  Below are the main functions used:

    ?BWInteracted   
    ?coefs  
    ?completeData  
    ?covStruct  
    ?ddfBW  
    ?ddfResidual  
    ?Ematrix    
    ?epsilon    
    ?hlTrace  
    ?Hmatrix  
    ?isReversible  
    ?LRT  
    ?MLtoREML  
    ?nestedCovariance  
    ?pillaiTrace  
    ?REMLtoML  
    ?royGR  
    ?sigStars  
    ?sphericityTests  
    ?termsType  
    ?uniTest  
    ?Vmatrix  
    ?waldF  
    ?wilksLambda  
    
The information on the datasets used can be accessed in the same way.

## Next Steps

This package only scratches the surface of marginal models and LMM.  It only addresses the first drawback up above.  Here are the next steps you can take to learn more about these types of models:

1) Fit different covariance structures other than unstructured and compound symmetry, and compare these covariance structures with likelihood ratio tests, AIC, and/or BIC with the `anova` function.

2) Learn about multiple comparisons and contrasts with the `emmeans` package.

3) Try out adding continuous time variables to your model.

4) Practice fitting competing models based on your research hypotheses and performing likelihood ratio tests for "reversible" and non-"reversible" models with the `anova` function.
