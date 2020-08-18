# A Spatial Bayesian Hierarchical model (BHM) for ground motion amplitudes
 
This model characterizes variations in earthquake ground motion parameters, accounting for event effects, station effects, event-station effects and unexplained effects. Previous models have assumed insignifiant ground motion variability at uniform site conditions over close distances, but observations have proven that is not the case. The motivation for this study is therefore to capture observed PGA variations and quantify to what extent the source and recording sites contribute to the overall variation in ground motions over relatively small distances on the same (lava-rock) site condition.

## The model

Detailed information about the project and the model can be found in the project's [research article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497). In this section we provide a summary of the model description and relevant parameters.

The proposed model: (Equation 1)

<img src="https://render.githubusercontent.com/render/math?math=Y_{es}=\mu_{es}(M_e,R_{es},D_e)%2B\delta+B_e%2B\delta+S_s%2B\delta+R_{es}"> 
<img src="https://render.githubusercontent.com/render/math?math=e=1,...,N,+s=1,...,Q">



with:

* <img src="https://render.githubusercontent.com/render/math?math=Y_{es}\text{: log10-PGA}">
* <img src="https://render.githubusercontent.com/render/math?math=\mu_{es}\text{: Median ground motion}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_e\text{: Effect from event e}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+B_e\text{: Effect from station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+S_s\text{: Effect from station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+WS_{es}\text{: Spatially correlated event-staion effect from event e and station s}">
* <img src="https://render.githubusercontent.com/render/math?math=\delta+R_{es}\text{: Effects that are unexplained or not accounted for}">


GMM for the median ground motion: (Equation 2)

<img src="https://render.githubusercontent.com/render/math?math=\mu_{es}=\beta_1%2B\beta_2+M_e%2B\beta_3+\text{log}_{10}(R_{es})%2B\beta_4+D_e">

with:

* <img src="https://render.githubusercontent.com/render/math?math=M_{e}\text{: local magnitude of the  }+e\text{th earthquake}">
* <img src="https://render.githubusercontent.com/render/math?math=\R_{es}\text{: Hypocentral distance from the  }+e\text{th event to the }s\text{th station}">
* <img src="https://render.githubusercontent.com/render/math?math=\D_e\text{: Depth of the origin of the  }+e\text{th earthquake}">

See [the article](https://onlinelibrary.wiley.com/doi/epdf/10.1002/env.2497) to see covariance matrix calculations and a further description of the model, parameters and the MCMC algorithm.


## The data

When setting up the model (see model setup file) you will receive a zip file containing 3 repositories and 6 Matlab files.

The repositories:

* Data: Containing all the input parameters
* Figs: Containing examples of how the figures will look like
* mat: Containing saved model outputs


The Matlab files:

* BHM_mainbody.m: The only file that you will run
* BHM_lnpost.m
* gpar.m
* sr_figr.m
* sr_hessian.m
* xcorr.m



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
* sim_latent: contains values for the 74 unknown latent parameters, for each iteration in each Markov chains



## Program setup

### Download Matlab

First of all, make sure that:
* Matlab version ___ or newer is installed on your computer  (This project was tested and created on Matlab version ____)
* Matlab's Parallell computing toolbox is installed (If installing Matlab for the first time, this can be installed simultaneously)

To go straight to the Matlab downloading site, press [here](https://nl.mathworks.com/downloads/).


### Download files
Second of all, download the whole project repository:
* Press the green "Code" button in the right top corner of the repository page
* Press Download ZIP

Extract the ZIP file to a folder at your preferred location.


### Insert PATH

So that Matlab can find the files on your computer, follow these steps:
* Open the recently downloaded BHM_mainbody.m
* In line 38 you can see the variable MAIN = 'C:/Users/...';
* Replace this line with the path to the project on your computer, i.e. MAIN='C:/Users/name_of_user/Documents/BHM_earthquake_gm'.

### Run

Now you are free to run the whole file, or run it section by section. If that is desired, click the section you want to run and press Ctrl+enter.


