# quantdr 1.2.1

* This is a resubmission.  The package was removed from CRAN due to dependency issues with MTS package.  
However, the dependence was removed when the 1.2.1 version was submitted.  The communication was not
clear and, although the dependency was resolved, the package was removed from CRAN.  

# quantdr 1.2.1

## Minor changes

* Remove the dependency with the MTS package

# quantdr 1.2.0

## Major changes

* Add the ValAR.R function for estimating the one-step ahead Value-at-Risk
* Adjust bandwidth selection for cqs.R and llqr.R to accommodate for Value-at-Risk calculation

## Minor changes

* Update the vignette to include explanations for the ValAR.R function
* Choose one bandwidth in case llqrcv.R outputs more (if there is a tie)
* Add the PerformanceAnalytics under Suggests in DESCRIPTION to accommodate for the example in the ValAR.R function
* Add tests for the ValAR.R function

# quantdr 1.1.0

## Major changes

* Fix error in llqr.R function when x0 is provided
* Replace the lm function in cqs.R to save computational time

## Minor changes

* Fix typos in vignette and update the plots using the ggplot2 package
* Update all the plots in README using the ggplot2 package
* Add the ggplot2 under Suggests in DESCRIPTION

# quantdr 1.0.0

* First release


