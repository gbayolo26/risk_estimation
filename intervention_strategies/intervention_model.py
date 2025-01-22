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
import pandas as pd


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
                            test_symptoms_type = 0, #0: to detect p proportion individuals with symptomps, 1:to detect individuals following a geometric distribution with parameter p
                            seed = 2, 
                            seed_open_ABM = 1,
                            days_of_quarantine = 100,
                            output_dir = "./simulation_results/",
                            save_every_iter = 10,
                            time_infection = 'estimation', #'constant' or True
                            quarantined_random_interactions = 0,
                            quarantine_adoption_fraction = 1,
                            tp_rate = 1,
                            tn_rate = 1):
            
    params = Parameters(input_param_file="../tests/data/baseline_parameters.csv",
                        param_line_number = 1,
                        output_file_dir="../data_test",
                        input_households="../tests/data/baseline_household_demographics.csv")
    
    #to set parameters in the open ABM
    params.set_param( "n_total", n_total)
    params.set_param( "end_time", t_end)
    params.set_param( "n_seed_infection", initial_inf)
    params.set_param( "rng_seed", seed_open_ABM)
    params.set_param( "quarantined_daily_interactions", quarantined_random_interactions)#1: quarantine only at work, 0: work and random
                           
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
    noise_quarantine = np.random.random(n_total) if quarantine_adoption_fraction < 1 else None
            
    # initialize the ranking method   
    ranked_method.init(n_total, initial_inf, t_0, test_availables, quarantine_household, test_household, p_SS, p_SM, seed, seed_open_ABM, day_since_infected_sympthoms, 
                       output_dir, age_group, house_no, time_infection, quarantined_random_interactions, test_symptoms_type, quarantine_adoption_fraction, tp_rate, tn_rate)
    
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
                
        if (t == t_0) & (seed !=2) & (test_symptoms_type == 0):            
            noise_SM = np.random.random(n_total) if p_SM < 1 else None
            noise_SS = np.random.random(n_total) if p_SS < 1 else None    
                    
        # to test individuals with severe symptomps   
        individuals_False_N_SS = []
        if test_symptoms_type == 0:               
                                        
             if p_SS < 1:  
                
               individuals_True_P_SS = individuals[(status == 4) & (noise_SS < p_SS)]
               
             else :
               individuals_True_P_SS = individuals[status == 4]   
                      
        if test_symptoms_type == 1:               
            
             individuals_SS = individuals[status == 4]
             
             if len(individuals_SS) > 0:
                                       
                 if p_SS < 1: 
                     noise = np.random.random(len(individuals_SS))
                     individuals_True_P_SS = individuals_SS[noise < p_SS]
                                     
                 else :
                     individuals_True_P_SS =  individuals_SS
                     individuals_False_N_SS = []
                     
             else :                 
                 individuals_True_P_SS =  []
                 
        #to test on individuals with mild symptomps   
        if test_symptoms_type == 0:               
            
             if tp_rate < 1: 
                 #tp_rate_mild = tp_rate + 0.05
                 tp_rate_mild = 1
                 individuals_True_P_SM = individuals[(status == 5) & (noise_SM <= p_SM*tp_rate_mild)]
                 individuals_False_N_SM = individuals[(status == 5) & (noise_SM > p_SM*tp_rate_mild)& (noise_SM < p_SM)]
            
             elif (p_SM < 1) & (tp_rate == 1):             
                 individuals_True_P_SM = individuals[(status == 5) & (noise_SM < p_SM)]
                 individuals_False_N_SM = []
                 
             else :
                 individuals_True_P_SM = individuals[status == 5]   
                 individuals_False_N_SM = []
        
        if test_symptoms_type == 1:               
            
             individuals_SM = individuals[status == 5]
                           
             if len(individuals_SM) > 0:
             
                 if tp_rate < 1:       
                     #tp_rate_mild = tp_rate + 0.05
                     tp_rate_mild = 1
                     noise = np.random.random(len(individuals_SM))
                     individuals_True_P_SM = individuals_SM[noise <= p_SM*tp_rate_mild]
                     individuals_False_N_SM = individuals_SM[(noise > p_SM*tp_rate_mild) & (noise < p_SM)]
            
                 elif (p_SM < 1) & (tp_rate == 1): 
                     noise = np.random.random(len(individuals_SM))
                     individuals_True_P_SM = individuals_SM[noise < p_SM]
                     individuals_False_N_SM = []
                 
                 else :
                     individuals_True_P_SM =  individuals_SM
                     individuals_False_N_SM = []
                     
                                 
             else :
                 
                 individuals_True_P_SM = []
                 individuals_False_N_SM = []
        
        individuals_True_P_SS =  np.setdiff1d(individuals_True_P_SS, all_TP)  
        
        individuals_True_P_SM =  np.setdiff1d(individuals_True_P_SM, all_TP)  
        
        #to upgrade the data test file
        ranked_method.test_observations_symptomps_upgrade(individuals_True_P_SS, individuals_False_N_SS, individuals_True_P_SM, individuals_False_N_SM, t)
        
        #number of tests at t
        test_availables_today = test_availables*(t >= t_0)   
        
        #to obtain a vector of ranked individuals        
        ranked_individuals = ranked_method.rank(t, daily_contacts)
       
        def f_test(individuals_to_test):
            
            s = status[individuals_to_test]            
            Pos = individuals_to_test[(s !=0) & (s < 8)]
            Neg = individuals_to_test[(s ==0) | (s >= 8)]
            
            if tp_rate < 1:
               noise = np.random.random(len(Pos))
               True_P = Pos[noise <= tp_rate] 
               False_N = Pos[noise > tp_rate] 
               
            else:
                True_P = Pos
                False_N = []    
                
            if tn_rate < 1:
               noise = np.random.random(len(Neg))
               True_N =  Neg[noise <= tn_rate]    
               False_P =  Neg[noise > tn_rate]   
               
            else: 
               True_N = Neg
               False_P = []    
                    
            Test_P = list(True_P) + list(False_P)
            Test_N = list(True_N) + list(False_N)
            
            return Test_P, Test_N, True_P, True_N        
       
        if len(ranked_individuals) > 0:
            
            #to select and test the individuals with the highest risk
            individuals_to_test_risk = np.array(ranked_individuals[ 0 : test_availables_today])
            TP_risk, TN_risk, True_P_risk, True_N_risk = f_test(individuals_to_test_risk)
            
            #to upgrade the data test file    
            ranked_method.test_observations_risk_upgrade(TP_risk, TN_risk, status, t)  
                    
        else:
            
            TP_risk = []
            TN_risk = []
            True_P_risk = []
        
        #individuals tested positive at t    
        #print(TP_risk)
        TP = list(TP_risk) + list( individuals_True_P_SM ) +  list(individuals_True_P_SS)  
        True_P = list(True_P_risk) + list(individuals_True_P_SM ) +  list(individuals_True_P_SS)  
                
        if quarantine_household == True:
           
           houses_infected = house_no[TP]
           idx = np.in1d(house_no, houses_infected)
           individuals_household = np.setdiff1d(individuals[idx], all_TP)   
           
           if test_household == True: 
               
              TP_house, TN_house, True_P_house, True_N_house = f_test(individuals_household)              
              ranked_method.test_observations_house_upgrade(TP_house, TN_house, status, t)    
              TP += list(TP_house)
              to_quarantine = TP
                               
           else:
               
              to_quarantine = TP + list(individuals_household)
              True_P_house = []
              
        else:    
            
           to_quarantine = TP
           True_P_house = []
        
        #to quarantine individuals at t
        if (quarantine_adoption_fraction < 1) & (len(to_quarantine) >=1):
            
            noise = noise_quarantine[to_quarantine]
            to_quarantine = np.array(to_quarantine)
            to_quarantine = list(to_quarantine[noise < quarantine_adoption_fraction])
        
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
        timeseries_sim['detected_R'].append(len(True_P_risk))
        timeseries_sim['detected_SS'].append(len(individuals_True_P_SS))     
        timeseries_sim['detected_SM'].append(len(individuals_True_P_SM))            
        timeseries_sim['detected_H'].append(len(True_P_house))        
        timeseries_sim['detected_active'].append(detected_active)
        timeseries_sim['n_quarantine'].append(len(to_quarantine)) 
        timeseries_sim['n_tests'].append(len(True_P))
        
        #print('time:', t, 'detected by risk:',len(TP_risk)) 
        
        #to save the results
        if t % save_every_iter == 0:            
            ranked_method.save_results(timeseries_sim, t)

    df_timeseries_sim = ranked_method.save_results(timeseries_sim, t_end)   
    
    sim.env.model.write_individual_file()
    df_indiv = pd.read_csv( "/home/gabriela/Desktop/OpenABM-Covid19-master(3)/OpenABM-Covid19-master/data_test/individual_file_Run1.csv", comment="#", sep=",", skipinitialspace=True )    
    ranked_method.save_individuals_results(df_indiv)   
    
    print('Total infected at time', t_end , 'for', ranked_method.name(), ':', n_total - Sucep)
    
    return df_timeseries_sim