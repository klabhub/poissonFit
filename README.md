# poissyFit

This Matlab code estimates the Poisson firing rate that underlies measured fluorescence in calcium imaging experiments.
 
 The theory behind this is described here:
  Ganmor, E., Krumin, M., Rossi, L. F., Carandini, M. & Simoncelli, E. P.
  Direct Estimation of Firing Rates from Calcium Imaging Data. 1–34 (2016).  http://arxiv.org/abs/1601.00364
  
  The key assumptions are:

  1. Spike counts are distributed according to a Poisson distribution with an underlying rate (lambda).
  2. Each spike produces a fixed amount of extra fluorescence (fPerSpike) and this decays exponentially over time (tau).
  3. The fluorescence is corrupted with additive Gaussian noise.

If these assumptions are met, this code will estimate a tuning curve (rate as a function of a set of conditions in an experiment)  better than the standard sequential approach in whcih  spike deconvolution is used to estimate spike times/counts, and then a tuning curve is estimated from these spike count estimates.  
  
  My implementation borrows heavily from the code provided at  https://github.com/eladganmor/Imaging_Analysis. I primarily added code comments, added optimization parameters, bootstrapping.
  
  
  Bart Krekelberg  - May 2023.
  
  
  ## Installation
  Clone this repository to your machine (e.g. to ```c:\github\poissonFit```), and add the folder to your Matlab path
  ```addpath('c:\github\poissonFit')```. 
  Then check out the demos folder:
  1. demos/simple.m : poisson rate estimation
  2. demos/parametric.m : parametric tuning curve estimation (von Mises)

  
  
  

