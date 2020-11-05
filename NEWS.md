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
