# Diffusion_Limited_Cancer_Growth_Model

Here we develop a mechanistic mathematical model to describe the growth and treatment of an avascular tumour for cancer drug development. We list both scripts for the non-spatial and spatial drug models. We also list the mixed effects scripts for the model without drug input. All the workflow is in Matlab.


This repository contains the following scripts:

Diffusion_Limited_Time_Dependent_Drug_Fitting_Function.m: This function solves the model with non-spatial drug input. The required inputs are parameters (a vector of 5 elements) and the timeframe of the experiment. The user can manually change the dose timings (currently set to dosing on days 1,8 and 15).

Two_Phase_Spatial_Drug_Function.m: This function solves the diffusion-limited model with spatial drug input. Includes plotting of drug concentration heat maps as well as the model solution.

Diffusion_Limited_Growth_Model_NLME.m: MatLab script to apply nonlinear mixed effects modelling to optimise parameters of the diffusion limioted model without drug input. This uses randomly generated mock data from the function Mock_Growth_Data_Generator.m described below. Uses the stochastic approximation to expectation maximisation 'nlmefitsa'.

Mock_Growth_Data_Generator.m: Randomly generates mock tumour growth data. This is done by drawing from a uniform distribution (from 0-1) and sorting into ascending order. Data for volume, time, subject and group are arranged into format for nlme fitting. 
			 
If using any parts of this code please cite

```
A. Nasim, J. Yates, G. Derks, C.M. A spatially resolved mechanistic growth law for cancer drug development predicting tumour growing fractions,
 Cancer Research Communications, 2022.
```
