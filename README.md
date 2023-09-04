# risk_estimation
## Test allocation based on risk of infection from first and second order contact tracing for Covid-19

This code contains four interventions methods that combines testing, contact tracing, and quarantine to mitigate the propagation of Covid-19,

 - First degree contact tracing
* Second degree contact tracing
+ Usual contact tracing
- Random selection

For more information see paper (to put link here). 
  
 We evaluate interventions based on our risk ranking through simulations of a fairly realistic agent-based model calibrated for COVID-19 epidemic outbreak ([the Oxford OpenABM model](https://github.com/BDI-pathogens/OpenABM-Covid19)). 
 This agent based model simulates the spread of the COVID-19 disease in a population with demographic structure based upon UK census data.  It has several implementations advantages, such as its level of detail and a very fast running time, even for large populations. 
 
 It permit to test containment strategies, in fact, contact tracing techniques such as the proposed can be easily integrated using the following version: [https://github.com/aleingrosso/OpenABM-Covid19](https://github.com/aleingrosso/OpenABM-Covid19). 

### Install
The code requires to install [the OpenABM-Covid19](https://github.com/aleingrosso/OpenABM-Covid19)

### Usage
The [test allocation notebook](https://github.com/gbayolo26/intervention_strategies/blob/main/Test_allocation.ipynb) contains an example of how to run the model
