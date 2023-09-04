#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 12:57:52 2023

@author: gabriela
"""

import numpy as np
import pandas as pd

def f_first_contacts(recent_contacts, data_test, t , gamma ):
    """
    This function compute the first degree interactions 
    
    """    
    all_TP = data_test["ID"][data_test["result"] == 1]               
    index_cases = data_test["ID"][(data_test['result'] == 1) & (data_test['time'] < t ) & (data_test['time'] >= t - gamma  )]
    
    df_FI = pd.DataFrame(recent_contacts, columns = ['ID', 'ID_2', 'time', 'type'])      
    df_FI  = df_FI[df_FI.ID.isin(index_cases ) & ~df_FI.ID_2.isin(all_TP) ] 
        
    return  df_FI


def f_number_of_contacts(df_FI, data_test):
    """
    This function compute the number of risky contacts by individual 
    
    """
    n_total = len(data_test['ID'])
    
    if len(df_FI)>0:
       #to compute the risk (number of contacts)
        df = df_FI[['ID', 'ID_2']].groupby(['ID_2'], as_index=False).count()
        df.columns = ['ID', 'Total_Risk']   
    
        #to put zero in the nan values
        df_Risk = pd.merge(pd.DataFrame({'ID': range(n_total)})  , df, how="outer")
        df_Risk['Total_Risk'] = df_Risk['Total_Risk'].fillna(0)
    else :
        df_Risk = pd.DataFrame({'ID': range(n_total), 'Total_Risk': 0})
        
    all_TP = data_test["ID"][data_test["result"] == 1]
    df_Risk.loc[df_Risk.ID.isin(all_TP), 'Total_Risk'] = -1  
    
    return df_Risk 
        

class Contact_Tracing:
        
    def __init__(self, gamma ):   
        self.description = "class for contact tracing ranking."        
        self.gamma = gamma  
            
    def init(self, n_total, initial_inf, t_0, test_availables, quarantine_household, test_household, p_SS, p_SM, seed, seed_open_ABM,
             day_since_infected_sympthoms, output_dir, age_group, house_no, real_time_infection, quarantined_random_interactions):
         
        self.n_total = n_total
        self.initial_inf = initial_inf
        self.t_0 = t_0
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
                
        # initialize first degree interactions 
        self.first_contacts = pd.DataFrame(columns = ['ID', 'ID_2', 'time', 'type']) 
        
        # initialize the recent interactions (interactions in the last gamma days)   
        self.recent_contacts = []
        
        #initialize the data test information     
        data_test = {}
        data_test["ID"] = np.arange(n_total)
        data_test["result"] = np.zeros(n_total) - 1
        data_test["time"] = np.zeros(n_total) - 1
        self.data_test = data_test
        
        np.random.seed(1)
        self.noisse = np.random.random(n_total)
           
        return True   
    
    def test_observations_symptomps_upgrade(self, individuals_SS, individuals_SM, t):
        
        self.data_test['time'][individuals_SS] = t            
        self.data_test['result'][individuals_SS] = 1            
        
        self.data_test['time'][individuals_SM] = t            
        self.data_test['result'][individuals_SM] = 1       
                
        return True
    
    def test_observations_house_upgrade(self, TP_house, TN_house, status, t):
         
        self.data_test['time'][TP_house] = t            
        self.data_test['result'][TP_house] = 1            
        
        return True
    
    def test_observations_risk_upgrade(self, TP_risk, TN_risk, status, t ):
        
        self.data_test['time'][TP_risk] = t            
        self.data_test['result'][TP_risk] = 1            
        
        return True
    
    def removals_upgrade(self,  status, t ):
        
        all_TP = np.array(self.data_test["ID"][self.data_test["result"] == 1])
        individuals_removed = list(all_TP[status[all_TP] >= 8])
        
        return len(individuals_removed)
   
    def rank(self, t, daily_contacts):
        '''
        to rank individuals given the number of contacts with index cases
        return: list of ranked individuals
        '''
        
        self.recent_contacts = [(col[0], col[1], col[2], col[3]) for col in self.recent_contacts if col[2] != t - self.gamma - 1]     
        self.recent_contacts.extend(daily_contacts)
           
        if t >= self.t_0:
             
            #to obtain the first order interactions 
            self.first_contacts = f_first_contacts(self.recent_contacts, self.data_test, t, self.gamma )
            
            #to compute the risk (number of contacts with index cases)
            df_Risk = f_number_of_contacts(self.first_contacts, self.data_test)                         
            df_Risk['noisse'] = np.random.random(1) - self.noisse
            df_Risk.loc[df_Risk['noisse'] <0, 'noisse'] = 2 + df_Risk.loc[df_Risk['noisse']<0, 'noisse']
            data_to_test = df_Risk[df_Risk['Total_Risk'] >= 0]                
            idx = np.lexsort((data_to_test['noisse'], data_to_test['Total_Risk']))
            
            #ranked individuals
            individuals_rank = list(data_to_test.iloc[idx]['ID'])
            individuals_rank.reverse()
                        
        else:
            
            individuals_rank = []    
            df_Risk = []
            
        return  individuals_rank
    
    def get_info(self):
        
        return self.data_test, self.first_contacts 
        
    def save_results(self, timeseries, t_end):        
                          
        name_file_res = 'timeseries_CT/timeseries_CT_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_0)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_qri_'+str(self.quarantined_random_interactions)
        name_data_test = 'timeseries_CT/data_test_CT_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_0)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_qri_'+str(self.quarantined_random_interactions)
        
        df_timeseries = pd.DataFrame(timeseries)                    
        df_timeseries['det_prop'] = (df_timeseries['detected_R'] + df_timeseries['detected_H'])/(df_timeseries['active'] - df_timeseries['detected_active'] + df_timeseries['detected_R']+ df_timeseries['detected_H'])
    
        if len(df_timeseries) < t_end:
    
            df_timeseries = pd.merge(pd.DataFrame({'time': range(t_end)}), df_timeseries, how="outer")
        
        df_timeseries.to_csv(self.output_dir + name_file_res+".csv")     
        dt = pd.DataFrame(self.data_test)
        dt.to_csv(self.output_dir + name_data_test+".csv")        
        
        return df_timeseries
        