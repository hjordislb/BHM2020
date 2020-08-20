# A Spatial Bayesian Hierarchical model (BHM) for ground motion amplitudes
 
This model characterizes variations in earthquake ground motion parameters, accounting for event effects, station effects, event-station effects and unexplained effects. Previous models have assumed insignifiant ground motion variability at uniform site conditions over close distances, but observations have proven that is not the case. The motivation for this study is therefore to capture observed PGA variations and quantify to what extent the source and recording sites contribute to the overall variation in ground motions over relatively small distances on the same (lava-rock) site condition.

## The model

Detailed information about the project and the model can be found in the project's [research article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497). In this section we provide a short summary of the model description and relevant parameters.

The proposed model: (Equation 1)

<img src="https://render.githubusercontent.com/render/math?math=log_{10}(Y_{es})=log_{10}(\mu_{es})%2B\delta+B_e%2B\delta+S_s%2B\delta+R_{es}"> 
<img src="https://render.githubusercontent.com/render/math?math=e=1,...,N,+s=1,...,Q">


with:

* <img src="https://render.githubusercontent.com/render/math?math=Y_{es}\text{: PGA}">
* <img src="https://render.githubusercontent.com/render/math?math=\mu_{es}\text{: Median ground motion}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_e\text{: Effect from event e}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_e\text{: Effect from station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+S_s\text{: Effect from station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+WS_{es}\text{: Spatially correlated event-staion effect from event e and station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+R_{es}\text{: Effects that are unexplained or not accounted for}">

Instead of modelling the peak ground acceleration (PGA), you could replace it with other ground motion parameters instead, i.e. peak ground velocity (PGV).

GMM for the median ground motion as a function of local magnitude, hypocentral distance and depth of the origin: (Equation 2)

<img src="https://render.githubusercontent.com/render/math?math=log_{10}(\mu_{es})=\beta_1%2B\beta_2+M_e%2B\beta_3+\text{log}_{10}(R_{es})%2B\beta_4+D_e">

with:

* <img src="https://render.githubusercontent.com/render/math?math=M_{e}\text{: local magnitude of the  }+e\text{th earthquake}">
* <img src="https://render.githubusercontent.com/render/math?math=\R_{es}\text{: Hypocentral distance from the  }+e\text{th event to the }s\text{th station}">
* <img src="https://render.githubusercontent.com/render/math?math=\D_e\text{: Depth of the origin of the  }+e\text{th earthquake}">

See [the article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497) to see covariance matrix calculations and a further description of the model, parameters and the MCMC algorithm.


## The data

When setting up the model (see model setup file) you will receive a zip file containing 3 repositories and 6 Matlab files. This data is provided in order for you to be able to test the model. It is not the original data from the research article.

The repositories:

* Data: Containing all the input parameters
* Figs: Containing examples of how the figures will look like
* mat: Containing saved model outputs


The Matlab files: (Not necessary to know their functionalities in order to run the model)

* BHM_mainbody.m: The main file - the only file that is manually run
* BHM_lnpost.m: This function is equivalent to the hyperparameter posterior function
* gpar.m: Computes the Gelman-Rubin statistic, the estimated effective sample size, the lag 1 autocorrelatin and the acceptance ratio.
* my_hessian.m: Computes the hessian matrix from the mode
* xcorr.m: Computes cross-correlation function estimates



### Input parameters

The input parameters for the model can be found in the repository "Data" 

* PGA_Obs: Ground motion amplitudes ( <img src="https://render.githubusercontent.com/render/math?math=Y_{es}"> observations)
* dij: Inter-station distance (used in covariance matrix calculations)
* M: Local magnitude (<img src="https://render.githubusercontent.com/render/math?math=M_{e}"> in model)
* Depth: Depth (<img src="https://render.githubusercontent.com/render/math?math=D_{e}"> in model)
* R_HYP: Hypocentral distance (<img src="https://render.githubusercontent.com/render/math?math=R_{es}"> in model)

Not used in our analysis:
* BAz: Back-azimuth i.e. direction
* R_EPI: Epicentral distance

### Output parameters

Since our model is a hierarchical model, the output parameters are divided into two parts: hyperparameters and latent parameters.

The hyperparameters are the following:
* <img src="https://render.githubusercontent.com/render/math?math=\tau">: inner-event std
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{S2S}">: interstation std
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{R}">: joint std in error term (unexplained events)
* <img src="https://render.githubusercontent.com/render/math?math=\phi_{SS}">: event-station std
* <img src="https://render.githubusercontent.com/render/math?math=\Delta_{SS}">: range parameter (within events)
* (<img src="https://render.githubusercontent.com/render/math?math=\Delta_{S2S}">: range parameter (station-to-station) - Fixed at 0.06 km) 

The latent parameters are:
* <img src="https://render.githubusercontent.com/render/math?math=\beta_1,...,\beta_4">, parameters from Equation 2
* <img src="https://render.githubusercontent.com/render/math?math=\delta+S_1,...,\delta+S_10">, parameters from Equation 1
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_1,...,\delta+B_60">, parameters from Equation 1


The output parameters for the model can be found in the repository "mat"

* sim_hyper: contains values for the 6 unknown hyperparameters, for each iteration in each Markov chain
* sim_latent: contains values for the 74 unknown latent parameters, for each iteration in each Markov chain.



## Program setup

### Install/Update Matlab

A few pointers
* This project was tested and created on Matlab version R2020a (Update 1). If there are problems with the project code, installing the newest version of Matlab might help. Testing on Matlab version R2019b was successful.
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
Extract the ZIP file to a folder at your preferred location.


### Insert PATH and run

So that Matlab can find the files on your computer, follow these steps:
* Open the recently downloaded BHM_mainbody.m
* Right at the start of the file (line 38) you can see the variable MAIN = 'C:/Users/...';
* Replace this line with the path to the project on your computer, i.e. MAIN='C:/Users/name_of_user/Documents/BHM_earthquake_gm'.
* Additionally, make sure that your "Current Folder" in Matlab is the project folder as well.

When this is done, you are free to run the whole project.


### Things to keep in mind

* Iterations: When the code is downloaded, the number of iterations in the Gibbs sampler will be set to NT=1.000 (in the "Gibbs sampler: Setup" section). This is done to make the testing of the code less time consuming. To acheive acceptable results, this parameter will need to be increased to 10.000-50.000.
* Mode sensitivity: The model is sensitive to the mode calculations. Sometimes a local mode is found instead of the global mode, which can alter the results. The figure to the left, down below, is an example of results using a global mode. On the figure to the right, the mode is local and does not give the correct results. If your plots look like the latter figure, find another mode and run the rest of the calculations again.

(Nota beta parameters og bara 1000 ítranir frekar)

![plots_hpar_50k](https://user-images.githubusercontent.com/39263646/90687955-378c5100-e25d-11ea-827e-5b572fef9d7b.jpg)

* Saving the mode: When a suitable mode has been found, you can comment out the mode finding section and fix that mode as the parameter "mode_theta" (First line in the "Gibbs sampler: Setup" section). This will save time later when the code is run again.

* Plot step size: In the "Convergence diagnostics" section of the code, this is the first line: II = 50:10:NT. The middle part determines the step size in the plots and will determine how time consuming the Convergence diagnostic section will be. When NT=10.000 or more, this could be changed to (i.e.) II = 50:100:NT to save time. 

* In the "Input" section of the main file the data is imported. The testing data is all .m files and therefore the code only loads .m files. Make sure to change the loading part so it is compatible with your file type. Loading of data of the following formats: .txt, .dat, or .csv, .xls, .xlsb, .xlsm, .xlsx, .xltm, .xltx, or .ods, can be found commented below the .m file loading part.



