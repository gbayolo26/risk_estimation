#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:13:00 2023

@author: gabriela
"""

import numpy as np

import covid19
from COVID19.model import Model, Parameters
import COVID19.simulation as simulation


def run_intervention_model( ranked_method, 
                            n_total = 50000, 
                            initial_inf = 10,
                            t_0 = 12,
                            t_end = 100, 
                            test_availables = 250, 
                            quarantine_household =  False, 
                            test_household = False,
                            p_SS = 1,
                            p_SM = 0.5,  
                            seed = 2, 
                            seed_open_ABM = 1,
                            days_of_quarantine = 100,
                            output_dir = "./simulation_results/",
                            save_every_iter = 10,
                            time_infection = 'estimation', #'constant' or True
                            quarantined_random_interactions = 0):
            
    params = Parameters(input_param_file="../tests/data/baseline_parameters.csv",
                        param_line_number = 1,
                        output_file_dir="../data_test",
                        input_households="../tests/data/baseline_household_demographics.csv")
    
    #to set parameters in the open ABM
    params.set_param( "n_total", n_total)
    params.set_param( "end_time", t_end)
    params.set_param( "n_seed_infection", initial_inf)
    params.set_param( "rng_seed", seed_open_ABM)
    params.set_param( "quarantined_daily_interactions", quarantined_random_interactions)#1 quarantine only at work, 0 work and random
                           
    #to iniatilize the open ABM model
    model = simulation.COVID19IBM(model = Model(params))
    sim   = simulation.Simulation(env = model, end_time = t_end, verbose=False)
      
    #to get information from the open ABM   
    day_since_infected_sympthoms = params.get_param("mean_time_to_symptoms") #get the mean value from the Open ABM    
    age_group = np.array(covid19.get_age(model.model.c_model))
    house_no = np.array(covid19.get_house(model.model.c_model))
        
    # initialize a fixed random stream before the intervention start    
    np.random.seed(2)
    noise_SM = np.random.random(n_total) if p_SM < 1 else None
    noise_SS = np.random.random(n_total) if p_SS < 1 else None
    
    # initialize the ranking method   
    ranked_method.init(n_total, initial_inf, t_0, test_availables, quarantine_household, test_household, p_SS, p_SM, seed, seed_open_ABM, day_since_infected_sympthoms, output_dir, age_group, house_no, time_infection, quarantined_random_interactions)
    
    # initialize the data    
    timeseries_sim = {}
    timeseries_sim['time'] = []    
    timeseries_sim['total_infected'] = []
    timeseries_sim['n_removed'] = []
    timeseries_sim['active'] = []
    timeseries_sim['detected_R'] = []    
    timeseries_sim['detected_SS'] = []
    timeseries_sim['detected_SM'] = []
    timeseries_sim['detected_H'] = []
    timeseries_sim['detected_active'] = []
    timeseries_sim['n_quarantine'] = []
    timeseries_sim['n_tests'] = []
    
    #variables to use after
    all_TP = [] 
    TP_house = []
    individuals = np.arange(n_total)
        
    for t in range(t_end):
        
        if (t == t_0) & (seed !=2):
            # initialize a fixed random stream when the intervention method starts   
            np.random.seed(seed)
            noise_SM = np.random.random(n_total) if p_SM < 1 else None
            noise_SS = np.random.random(n_total) if p_SS < 1 else None
                    
        #to run the Open ABM model one day
        sim.steps(1)  
        
        #to extract the individuals status
        status = np.array(covid19.get_state(model.model.c_model))                  
        Sucep = (status == 0).sum()
        Removed = (status >=8).sum()
        Infected = n_total - Sucep - Removed
               
        if Infected == 0:
            print("BREAK: no more infected individuals")            
            timeseries_sim['time'].append(t)
            timeseries_sim['total_infected'].append(n_total - Sucep)
            timeseries_sim['n_removed'].append(Removed)
            timeseries_sim['active'].append(Infected)
            timeseries_sim['detected_R'].append(0)
            timeseries_sim['detected_SS'].append(0)     
            timeseries_sim['detected_SM'].append(0)            
            timeseries_sim['detected_H'].append(0)        
            timeseries_sim['detected_active'].append(0)
            timeseries_sim['n_quarantine'].append(0) 
            timeseries_sim['n_tests'].append(0)
            
            break
        
        #to upgrade individuals removals (upgrade in the data test file the time of removal and infection)
        n_removals = ranked_method.removals_upgrade(status, t )
                
        #to extract daily contacts
        daily_contacts = covid19.get_contacts_daily(model.model.c_model, t)
                
        # to test individuals with severe symptomps                     
        if p_SS < 1:            
            individuals_SS = individuals[(status == 4) & (noise_SS < p_SS)]
        else :
            individuals_SS = individuals[status == 4]   
            
        individuals_SS =  np.setdiff1d(individuals_SS, all_TP)  
        
        #to test on individuals with mild symptomps        
        if p_SM < 1:            
            individuals_SM = individuals[(status == 5) & (noise_SM < p_SM)]
        else :
            individuals_SM = individuals[status == 5]   
            
        individuals_SM =  np.setdiff1d(individuals_SM, all_TP)  
        
        #to upgrade the data test file
        ranked_method.test_observations_symptomps_upgrade(individuals_SS, individuals_SM, t)
        
        #number of tests at t
        test_availables_today = test_availables*(t >= t_0)   
        
        #to obtain a vector of ranked individuals        
        ranked_individuals = ranked_method.rank(t, daily_contacts)
       
        def f_test(individuals_to_test):
            
            s = status[individuals_to_test]            
            TP = individuals_to_test[(s !=0) & (s < 8)]
            TN = individuals_to_test[(s ==0) | (s >= 8)]
                       
            return TP, TN        
       
        if len(ranked_individuals) > 0:
            
            #to select and test the individuals with the highest risk
            individuals_to_test_risk = np.array(ranked_individuals[ 0 : test_availables_today])
            TP_risk, TN_risk = f_test(individuals_to_test_risk)
            
            #to upgrade the data test file    
            ranked_method.test_observations_risk_upgrade(TP_risk, TN_risk, status, t)  
                    
        else:
            
            TP_risk = []
            TN_risk = []
        
        #individuals tested positive at t    
        TP = list(TP_risk) + list( individuals_SM ) +  list(individuals_SS)  
        
        if quarantine_household == True:
           
           houses_infected = house_no[TP]
           idx = np.in1d(house_no, houses_infected)
           individuals_household = np.setdiff1d(individuals[idx], all_TP)   
           
           if test_household == True: 
               
              TP_house, TN_house = f_test(individuals_household)              
              ranked_method.test_observations_house_upgrade(TP_house, TN_house, status, t)    
              TP += list(TP_house)
              to_quarantine = TP
                               
           else:
               
              to_quarantine = TP + list(individuals_household)
              
        else:    
            
           to_quarantine = TP
        
        #to quarantine individuals at t
        covid19.intervention_quarantine_list(model.model.c_model, to_quarantine, days_of_quarantine)
                 
        if time_infection == True:           
           ranked_method.real_time_infection_upgrade(TP)
           
        if time_infection == 'constant':
           ranked_method.constant_time_infection_upgrade(t, TP_risk)
           
        all_TP += TP
        
        detected_active = len(all_TP) - n_removals        
        timeseries_sim['time'].append(t)
        timeseries_sim['total_infected'].append(n_total - Sucep)
        timeseries_sim['n_removed'].append(Removed)
        timeseries_sim['active'].append(Infected)
        timeseries_sim['detected_R'].append(len(TP_risk))
        timeseries_sim['detected_SS'].append(len(individuals_SS))     
        timeseries_sim['detected_SM'].append(len(individuals_SM))            
        timeseries_sim['detected_H'].append(len(TP_house))        
        timeseries_sim['detected_active'].append(detected_active)
        timeseries_sim['n_quarantine'].append(len(to_quarantine)) 
        timeseries_sim['n_tests'].append(len(TP))
        
        #print('time:', t, 'detected by risk:',len(TP_risk)) 
        
        #to save the results
        if t % save_every_iter == 0:            
            ranked_method.save_results(timeseries_sim, t)

    df_timeseries_sim = ranked_method.save_results(timeseries_sim, t_end)   
    
    
    
    print('CT method for seed', seed)
    
    return df_timeseries_sim