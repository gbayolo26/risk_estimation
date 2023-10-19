#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 17:16:55 2023

@author: gabriela
"""
import numpy as np
import pandas as pd

class Random_Selection:
        
    def __init__(self):
        self.description = "class for random selection method."        
        
    def init(self, n_total, initial_inf, t_start, test_availables, quarantine_household, test_household, p_SS, p_SM, seed, seed_open_ABM, 
             day_since_infected_sympthoms, output_dir, age_group, house_no, real_time_infection, quarantined_random_interactions):
         
        self.n_total = n_total
        self.initial_inf = initial_inf
        self.t_start = t_start
        self.test_availables = test_availables
        self.quarantine_household = quarantine_household
        self.test_household = test_household
        self.p_SS = p_SS
        self.p_SM = p_SM
        self.seed = seed
        self.seed_open_ABM = seed_open_ABM
        self.day_since_infected_sympthoms = day_since_infected_sympthoms
        self.output_dir = output_dir
        self.quarantined_random_interactions = quarantined_random_interactions        
        data_test = {}
        data_test["ID"] = np.arange(n_total)
        data_test["result"] = np.zeros(n_total) - 1
        data_test["time"] = np.zeros(n_total) - 1        
        self.data_test = data_test        
        np.random.seed(1)
        
        return True   
    
    def test_observations_symptomps_upgrade(self, individuals_SS, individuals_SM, t):
        
        self.data_test['time'][individuals_SS] = t            
        self.data_test['result'][individuals_SS] = 1            
        
        self.data_test['time'][individuals_SM] = t            
        self.data_test['result'][individuals_SM] = 1       
                
        return True
    
    def test_observations_risk_upgrade(self, TP_risk, TN_risk, status, t ):
        
        self.data_test['time'][TP_risk] = t            
        self.data_test['result'][TP_risk] = 1            
        
        return True
    
    def test_observations_house_upgrade(self, TP_house, TN_house, status, t):
         
        self.data_test['time'][TP_house] = t            
        self.data_test['result'][TP_house] = 1            
        
        return True
    
    def removals_upgrade(self,  status, t ):
        
        all_TP = np.array(self.data_test["ID"][self.data_test["result"] == 1])
        individuals_removed = list(all_TP[status[all_TP] >= 8])
        
        return len(individuals_removed)
   
    def rank(self, t, daily_contacts):
        '''
        to rank individuals randomly
        return: list of ranked individuals
        '''
        
        if t >= self.t_start:
            
            individuals_rank = np.random.choice( list(self.data_test['ID'][self.data_test['result'] != 1]), self.test_availables, replace = False)
                                    
        else:
            
            individuals_rank = []    
            
        return  individuals_rank
    
    def name(self):
        
        return 'RS'
    
    def save_results(self, timeseries, t_end):        
                          
        name_file_res = 'timeseries_RS_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_qri_'+str(self.quarantined_random_interactions)   
        name_data_test = 'data_test_RS_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_qri_'+str(self.quarantined_random_interactions)   
         
        df_timeseries = pd.DataFrame(timeseries)
        
        if len(df_timeseries) < t_end:
    
            df_timeseries = pd.merge(pd.DataFrame({'time': range(t_end)}), df_timeseries, how="outer")
        
        df_timeseries.to_csv(self.output_dir + name_file_res+".csv")    
        dt = pd.DataFrame(self.data_test)
        dt.to_csv(self.output_dir + name_data_test+".csv")        
        
        return df_timeseries
        
    
    
        