## Changes in Version 2.1.0
  - Added a more robust solve method for armadillo matrices
  - Fixed issues with exported Bhattacharya and Rue algorithm 
  - Fixed issue with the usenames option for xcolor in the vignette

## Changes in Version 2.0.6
  - Fixed a bug in LPDS()
  - Fixed issues on Fedora and Debian 
  - Changed compiler flag -DARMA_DONT_PRINT_ERRORS to -DARMA_WARN_LEVEL=0

## Changes in Version 2.0.5
  - Modified shrinkTVP.ltx to work with new jss.cls

## Changes in Version 2.0.4
  - Added CITATION file
  - Modified DESCRIPTION as well as some documentation files

## Changes in Version 2.0.3
  - Fixed a bug in the sampling step for c_xi (Thanks to Tony Chernis for pointing it out)
  - Fixed a bug in the sampling step for a_tau (Thanks to Wenjie Zhao for pointing it out)
  - Added DG_MH_step to exported functions
  - Fixed some typos in the documentation
  - Added a newline to the print() function

## Changes in Version 2.0.2
  - Moved to stochvol 3.0
  - Fixed an issue with the exported version sample_TG_TVP function
  - Fixed an issue with the exported Bhattacharya algorithm 
  - define STRICT_R_HEADERS, include float.h, adjust some constants

## Changes in Version 2.0.1
  - Fixed an issue on Solaris

## Changes in Version 2.0
  - Added option to use triple gamma prior 
  - Added one-step update function updateTVP for modular use in other samplers
  - Added convenience methods (residuals.shrinkTVP, fitted.shrinkTVP, etc...)
  - Re-factored C++ code, improving on speed in some places
  - Added forecast_shrinkTVP (and associated plot.forecast_shrinkTVP method), which allows for forecasting
  - All C++ functions are now exported for use in other packages
  - simTVP now allows user to adjust parameters governing SV process
  - Re-factored plot.shrinkTVP to use layout() instead of par(mfrow = ...)
  - Bug fixes
  - Updated the vignette

## Changes in Version 1.1.1
  - Updated the vignette.

## Changes in Version 1.1.0
  - Overhauled plotting functions. The default is now shaded areas instead of lines and more customisation options are available.
  - plot.mcmc.tvp now also adds time series information to the x-axis labels, thanks to zoo::index.
  - Changed the way LPDS are calculated, moving them from the shrinkTVP function to a dedicated LPDS function.
  - Added the new eval_pred_dens function, which allows the user to evaluate the one-step ahead predctive density.
  - Objects returned from the shrinkTVP function are now of class shrinkTVP and not shrinkTVP_res. 

## Changes in Version 1.0.2:
  - Added a vignette
  - Bug fixes


## Changes in Version 1.0.1:
  - plot.mcmc.tvp now adds a horizontal line at 0
  - Added a nicer print method for shrinkTVP_res objects
  - Fixed some typos and errors in the documentation
