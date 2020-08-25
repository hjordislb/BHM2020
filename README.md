# A Spatial Bayesian Hierarchical model (BHM) for ground motion amplitudes
 
This model characterizes variations in earthquake ground motion parameters, accounting for event effects, station effects, event-station effects and unexplained effects. Previous models have assumed insignifiant ground motion variability at uniform site conditions over close distances, but observations have proven that is not the case. The motivation for this model is therefore to capture observed ground motion variations (measured in peak ground acceleration) and quantify to what extent the source and recording sites contribute to the overall variation in ground motions over relatively small distances on the same (lava-rock) site condition.

## The model

Detailed information about the project and the model can be found in the project's [research article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497). In this section we provide a short summary of the model description and relevant parameters.

The proposed model: (Equation 1)

<figure>
  <img src="https://render.githubusercontent.com/render/math?math=\text{log}+_{10}(Y_{es})=\text{log}_{10}(\mu_{es})%2B\delta+B_e%2B\delta+S_s%2B\delta+WS_{es}%2B\delta+R_{es}"> 
  <img src="https://render.githubusercontent.com/render/math?math=e=1,...,N,+s=1,...,Q">
</figure>


with:

* N: Maximum number of events
* Q: Maximum number of stations

and:

* <img src="https://render.githubusercontent.com/render/math?math=Y_{es}\text{: Peak ground acceleration (PGA) recorded for event e and station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\mu_{es}\text{: Median ground motion}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_e\text{: Effect from event e}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+S_s\text{: Effect from station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+WS_{es}\text{: Spatially correlated event-staion effect from event e and station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+R_{es}\text{: Effects that are unexplained or not accounted for}">


Instead of modelling the peak ground acceleration (PGA), you could replace it with other ground motion parameters instead, i.e. peak ground velocity (PGV).

To obtain the median ground motion, it is modelled as a function of local magnitude, hypocentral distance and depth of the origin: (Equation 2)

<img src="https://render.githubusercontent.com/render/math?math=\text{log}+_{10}(\mu_{es})=\beta_1%2B\beta_2+M_e%2B\beta_3+\text{log}_{10}(R_{es})%2B\beta_4+D_e">

with:

* <img src="https://render.githubusercontent.com/render/math?math=M_{e}\text{: local magnitude of event e}">
* <img src="https://render.githubusercontent.com/render/math?math=\R_{es}\text{: Hypocentral distance from event e to station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\D_e\text{: Depth of the origin of event e}">


See [the article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497) to see covariance matrix calculations and a further description of the model, parameters and the MCMC algorithm.


## The data

When setting up the model (see "Program setup" down below) you will receive a zip file containing 3 repositories and 6 Matlab files. This data is provided in order for you to be able to test the model. It is not the original data from the research article.

The repositories:

* Data: Contains all the input parameters
* Figs: Contains examples of how the figures should look like after 50.000 gibbs-loop iterations
* mat: Contains saved model parameter outputs (will be overwritten after every run)


The Matlab files: (Not necessary to know their functionalities in order to run the model)

* BHM_mainbody.m: The main file (Equations are numbered with the corresponding equations in the research article)
* BHM_lnpost.m: This function is equivalent to the log-posterior of the input hyperparameters
* gpar.m: Computes the Gelman-Rubin statistic, the estimated effective sample size, the lag 1 autocorrelatin and the acceptance ratio.
* my_hessian.m: Computes the hessian matrix from the input mode
* xcorr.m: Computes cross-correlation function estimates
* randnLimit: Generates random normal values within input interval


### Input parameters

The input parameters to test the model can be found in the repository "Data" 

* PGA_Obs: Ground motion amplitudes ( <img src="https://render.githubusercontent.com/render/math?math=Y_{es}"> observations)
* dij: Inter-station distance (used in covariance matrix calculations)
* M: Local magnitude (<img src="https://render.githubusercontent.com/render/math?math=M_{e}"> in model)
* Depth: Depth (<img src="https://render.githubusercontent.com/render/math?math=D_{e}"> in model)
* R_HYP: Hypocentral distance (<img src="https://render.githubusercontent.com/render/math?math=R_{es}"> in model)


### Output parameters

Since our model is a hierarchical model, the output parameters are divided into two parts: hyperparameters and latent parameters.

The hyperparameters are the following:
* <img src="https://render.githubusercontent.com/render/math?math=\tau">: inter-event variability
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{S2S}">: inter-station variability
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{R}">: unexplained variability
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{SS}">: event-station variability
* <img src="https://render.githubusercontent.com/render/math?math=\Delta_{SS}">: range parameter (within events)
* (<img src="https://render.githubusercontent.com/render/math?math=\Delta_{S2S}">: range parameter (station-to-station) - Fixed at 0.06 km) 

The latent parameters are:
* <img src="https://render.githubusercontent.com/render/math?math=\beta_1,...,\beta_4">: ground motion model coefficients
* <img src="https://render.githubusercontent.com/render/math?math=\delta+S_1,...,\delta+S_10">: site-to-site amplification terms (inter-station residuals) 
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_1,...,\delta+B_60">: event terms (inter-event residuals)

The "delta-S" parameters represent the ten stations where the ground motions are measured. The "delta-B" parameters represent the 60 events of ground motion following an earthquake (in this case, within a year of the original earthquake).

The output parameters will be saved in the repository "mat":

* sim_hyper: contains values for the 6 unknown hyperparameters, for each iteration in each Markov chain
* sim_latent: contains values for the 74 unknown latent parameters, for each iteration in each Markov chain.

The plots will be saved in the repository "Figs" in .m format. All parameters and plots that are saved after running the code will be overwritten in the next run.


## Program setup

### Install/Update Matlab

* This project was tested and created in Matlab version R2020a (Update 1). If there are problems with the project code, installing the newest version of Matlab might help. Testing on Matlab version R2019b was also successful.
* Matlab's Parallel Computing Toolbox is used in the project and therefore has to be installed before running the code (If installing Matlab for the first time, this can be installed simultaneously)

For those who have not installed Matlab yet, follow [these](https://nl.mathworks.com/help/install/) steps to download and run the installer.
For those who have Matlab installed but do not have the Parallel Computing Toolbox. Follow these steps:

* Open the Home tap in Matlab
* Press Add-Ons, then Get Add-Ons
* Search for Parallel Computing Toolbox
* Open and install it


### Download files
The next step is to download the whole project repository:
* Press the green "Code" button in the right top corner of the GitHub repository page
* Press Download ZIP
* Extract the ZIP file to a folder at your preferred location.


### Insert PATH and run

 Follow these steps to make sure that Matlab can find the files on your computer:
* Open the recently downloaded BHM_mainbody.m
* Right at the start of the file (line 40) you can see the variable MAIN = 'C:/Users/...';
* Replace this line with the path to the project on your computer, i.e. MAIN='C:/Users/name_of_user/Documents/BHM2020'.
* Additionally, make sure that your "Current Folder" in Matlab is the project folder as well.

When this is done, you are free to run the BHM_mainbody file.


### Things to keep in mind

* Iterations: When the code is downloaded, the number of iterations in the Gibbs sampler will be set to NT=1.000 (in the "Gibbs sampler: Setup" section). This is done to make the testing of the code less time consuming. To acheive acceptable results, this parameter will need to be increased to 10.000-50.000.

* Mode sensitivity: The model is sensitive to the mode calculations. Sometimes a local mode is found instead of the global mode, which can alter the results. Figure 1, down below, is an example of results (after 1000 iterations) using a global mode. In Figure 2 the mode is local and will not give the correct results. If your beta plots look like the latter figure when using NT=1000, i.e. includes a lot of straight lines, find another mode and run the rest of the calculations again.

<figure>
  <figcaption>Figure 1</figcaption>
  <img src="https://user-images.githubusercontent.com/39263646/90778821-a1077080-e2ec-11ea-8874-429a7318e08e.jpg" width="600" height="400" />
</figure>

<figure>
  <figcaption>Figure 2</figcaption>
  <img src="https://user-images.githubusercontent.com/39263646/90778901-bb414e80-e2ec-11ea-9006-c8973e9057e4.jpg" width="600" height="400" />
</figure>

* Saving the mode: When a suitable mode has been found, you can comment out the "Global mode optimization" section and fix that mode as the parameter "mode_theta" (Can be done in the "Fixing the mode" section). This will save time later when the code is run again.

* Sensitivity analysis: When testing the model on your own data (not the testing data) some parameters might need to be tweaked. In particular: 
 - NT (Number of iterations in Gibbs-loop)
 - NC (Number of chains)
 - DeltaS2S (The range parameter fixed at 0.6 km)
In the testing data, sensitivity analysis implied that the best results (within reasonable time limits) were found using: NT=10.000, NC=10 and DeltaS2S=0.6. However, these values are not necessarily suitable for other datasets.

* In the "Input" section of the main file the data is imported. The testing data can be found as .m files and .xlsx files but the code is programmed to load .m files. Make sure to change the loading part, if necessary, so it is compatible with your file type. Loading of data of the following formats: .txt, .dat, or .csv, .xls, .xlsb, .xlsm, .xlsx, .xltm, .xltx, or .ods, can be found commented below the .m file loading part.



