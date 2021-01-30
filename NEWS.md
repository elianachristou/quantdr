# quantdr 2.1.0

## Major changes

* Add the ValAR.R function for estimating the value-at-Risk
* Update the vignette to include explanations for the ValAR.R function

## Minor changes

* Choose one bandwidth in case llqrcv.R outputs more (if there is a tie)
* Adjust bandwidth selection for cqs.R and llqr.R to accomodate for Value-at-Risk calculation

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


