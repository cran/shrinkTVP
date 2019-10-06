## Changes in Version 1.1.0
  - Overhauled plotting functions. The default is now shaded areas instead of lines and more customisation options are available.
  - plot.mcmc.tvp now also adds time series information to the x-axis labels, thanks to zoo::index.
  - Changed the way LPDS are calculated, moving them from the shrinkTVP function to a dedicated LPDS function.
  - Added the new eval_pred_dens function, which allows the user to evaluate the one-step ahead predctive density.
  - Objects returned from the shrinkTVP function are now of class shrinkTVP and not shrinkTVP_res. 
  - Updated the vignette.

## Changes in Version 1.0.2:
  - Added a vignette
  - Bug fixes


## Changes in Version 1.0.1:
  - plot.mcmc.tvp now adds a horizontal line at 0
  - Added a nicer print method for shrinkTVP_res objects
  - Fixed some typos and errors in the documentation
