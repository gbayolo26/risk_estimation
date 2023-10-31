# risk_estimation
## Test allocation based on risk of infection from first and second order contact tracing for Covid-19

This repository contains four interventions methods to mitigate the propagation of Covid-19, 

 - First degree contact tracing  ([code here](https://github.com/gbayolo26/risk_estimation/blob/main/intervention_strategies/First_Degree_Contact_Tracing_Risk.py))
* Second degree contact tracing  ([code here](https://github.com/gbayolo26/risk_estimation/blob/main/intervention_strategies/Second_Degree_Contact_Tracing_Risk.py))
+ Usual contact tracing ([code here](https://github.com/gbayolo26/risk_estimation/blob/main/intervention_strategies/Contact_Tracing.py))
- Random selection ([code here](https://github.com/gbayolo26/risk_estimation/blob/main/intervention_strategies/Random_Selection.py))

For more information see paper (to put link here). 
  
 We evaluate the non pharmaceutical interventions based on our risk ranking through simulations of a fairly realistic agent-based model calibrated for COVID-19 epidemic outbreak ([the Oxford OpenABM model](https://github.com/BDI-pathogens/OpenABM-Covid19)). 
 This agent based model simulates the spread of the COVID-19 disease in a population with demographic structure based upon UK census data.  
 
 It permit to test containment strategies, in fact, contact tracing techniques such as the proposed can be easily integrated using the following version: [https://github.com/aleingrosso/OpenABM-Covid19](https://github.com/aleingrosso/OpenABM-Covid19). 

### Install
The code requires to install [the OpenABM-Covid19](https://github.com/aleingrosso/OpenABM-Covid19)

### Usage
The [test allocation notebook](https://github.com/gbayolo26/intervention_strategies/blob/main/Test_allocation.ipynb) contains an example of how to run the model

### Maintainers
Gabriela Bayolo Soler (gabriela.bayolo-soler@utc.fr)

<img src="![utc](https://github.com/gbayolo26/risk_estimation/assets/79975920/eedd8e6e-6cea-4327-bef3-54fe0256ff06)" width="160" height="50">

### Licence
[Apache-2.0 licence](https://github.com/gbayolo26/risk_estimation/blob/main/LICENSE)
