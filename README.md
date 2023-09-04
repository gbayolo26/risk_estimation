# intervention_strategies
## Test allocation based on risk of infection from first and second order contact tracing for Covid-19

Strategies such as testing, contact tracing, and quarantine have been proven to be essential mechanisms to mitigate the propagation of infectious diseases.  
However, when an epidemic spreads rapidly and/or the resources to contain it are limited (e.g. not enough daily available tests), to test and quarantine all the contacts of detected individuals is impracticable.
In this direction, we propose a method to compute the individual risk of infection over time, based on the partial observation of the epidemic spreading through the population contact network. 
 We define the risk of an individual as her/his probability of getting infected from any of the possible chains of transmission up to length two, originating from recently detected individuals.
 We provide explicit formulas and an efficient implementation of our method.
 Ranking individuals according to their risk of infection can be used as a decision tool to determine which individuals get tested, quarantined, or applied other preventive measures in priority. 
 We evaluate interventions based on our risk ranking through simulations of a fairly realistic agent-based model calibrated for COVID-19 epidemic outbreak (the Oxford OpenABM model). We consider different scenarios to study the role of key quantities such as the number of daily available tests, the contact tracing time-window, the transmission probability per contact (constant versus depending on multiple factors), and the age since infection (for varying infectiousness). 
 We find that, when there is a limited number of daily tests available, our method is capable to  mitigate the propagation more efficiently than random selection, than the usual contact tracing (ranking according to the number of contacts with detected individuals), and than some other approaches in the recent literature on the subject. We stress, in particular, the important role played by the estimation (or knowledge) of the time since infection of index cases, in order to trace accurately their most risky contacts.
