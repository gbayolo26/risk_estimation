#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 10:26:10 2023

@author: gabriela
"""
 
import numpy as np
import pandas as pd
from.utils import f_compute_transmission_probability, f_time_to_detection, f_compute_Sa

def f_first_interactions(recent_contacts, data_test, t , prob_function, Sa):
   
    """
    This function compute the 1°interactions,
         
    inputs:
        recent_contacts : The recent interactions from the Open ABM Model
        data_test : Data frame with test information
        t : day of intervention
        
    outputs:
        df_FI : The first degree intercations 
        
    """
    
    all_TP = data_test[(data_test['result'] == 1) & (data_test['time'] < t )]["ID"]
    
    infected_period = 20
    time_active_test = list(range(max(0, t - infected_period), t))
    index_cases = data_test[data_test.result.isin([1]) & data_test.pos_tim_inf.isin(time_active_test) & ~data_test.time.isin([t])]['ID'] 
        
    df_FI = pd.DataFrame(recent_contacts, columns = ['ID', 'ID_2', 'time', 'type'])     
    df_FI  = df_FI[df_FI.ID.isin(index_cases ) & ~df_FI.ID_2.isin(all_TP) ] 
    
    if len(df_FI) > 0:
        
        df_FI = df_FI.merge(data_test[['ID', 'pos_tim_inf', 'time_removed', 'observed_status']].rename(columns = { 'pos_tim_inf':'pos_tim_inf_1', 
                                                                                                        'time_removed':'time_removed_1',
                                                                                                        'observed_status':'observed_status_1'}))
        
        df_FI  = df_FI [df_FI['pos_tim_inf_1'] < df_FI ['time']] 
        df_FI  = df_FI [df_FI['time_removed_1'] > df_FI ['time']] 
       
        df_FI = df_FI.merge(data_test[['ID', 'time', 'age_group']].rename(columns = { 'ID': 'ID_2', 'time':'time_test_neg_2', 'age_group':'age_group_2'}))
            
        df_FI  = df_FI [df_FI['time_test_neg_2'] < df_FI ['time']] 
        
    if len(df_FI) > 0:           
                    
        df_FI['Prob'] = f_compute_transmission_probability(df_FI, prob_function, Sa)
        df_FI['CProb'] = 1 - df_FI['Prob']
        
    else :
        df_FI = pd.DataFrame(columns = ['ID', 'ID_2', 'time', 'type'])
       
    return  df_FI

def f_first_contacts_risk(df_FI, data_test):
    
    """
    This function compute the risk for first degree contacts
    
    """
    
    n_total = len(data_test['ID'])
    
    if len(df_FI)>0:
        
        #to compute the first degree contact tracing risk    
        df = df_FI.groupby(['ID_2'], as_index=False).agg({'CProb': 'prod'})
        df.columns = ['ID', 'Total_Risk'] 
        df['Total_Risk'] = 1 - df['Total_Risk']    
    
        #to put zero in the nan values
        df_Risk = pd.merge(pd.DataFrame({'ID': range(n_total)})  , df, how="outer")
        df_Risk['Total_Risk'] = df_Risk['Total_Risk'].fillna(0)
    else :
        df_Risk = pd.DataFrame({'ID': range(n_total), 'Total_Risk': 0})
    
    all_TP = data_test[data_test["result"] == 1]["ID"]    
    df_Risk.loc[df_Risk.ID.isin(all_TP), 'Total_Risk'] = -1  
    
    return df_Risk 
 

class First_Degree_Contact_Tracing_Risk:
        
    def __init__(self, gamma, prob_function):
        
        self.description = "class for first degree contact tracing method."        
        self.gamma = gamma 
        self.prob_function = prob_function  
        
    def init(self, n_total, initial_inf, t_start, test_availables, quarantine_household, test_household, p_SS, p_SM, seed, seed_open_ABM,
             day_since_infected_sympthoms, output_dir, age_group, house_no, time_infection, quarantined_random_interactions):
         
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
        self.time_infection = time_infection 
        self.quarantined_random_interactions = quarantined_random_interactions
        self.day_since_infected_sympthoms_SS = day_since_infected_sympthoms 
        self.day_since_infected_sympthoms_SM = day_since_infected_sympthoms 
        self.recent_contacts = []
        self.first_interactions = []
        self.individuals_removed = list()
        
        if self.prob_function == 'lambda':
           self.Sa = f_compute_Sa(age_group, house_no)
           
        else:
            self.Sa = []
        
        data_test = pd.DataFrame()
        data_test["ID"] = np.arange(n_total)
        data_test["result"] = np.zeros(n_total) - 1
        data_test["type"] = np.zeros(n_total) - 1 #0: detection by symptomps, 1: by risk, 2: by home, -1: not detected
        data_test["time"] = np.zeros(n_total) - 2
        data_test["pos_tim_inf"] = np.zeros(n_total) - 1
        data_test["observed_status"] = np.zeros(n_total) 
        data_test["age_group"] = age_group
        data_test["house_no"] = house_no
        data_test["time_removed"] = np.zeros(n_total) + np.inf
        data_test["time_infected"] = np.zeros(n_total) - 1
        data_test["type"] = np.zeros(n_total) - 1
        self.data_test = data_test
        
        np.random.seed(1)
        self.noisse = np.random.random(n_total)
        
        return True
        
    def test_observations_symptomps_upgrade(self, individuals_SS, individuals_SM, t):
        
        idx_test_positive_SS = self.data_test.ID.isin(individuals_SS)        
        self.data_test.loc[idx_test_positive_SS,'pos_tim_inf'] = [max(1 + self.data_test["time"][i], t - self.day_since_infected_sympthoms_SS) for i in individuals_SS]
        self.data_test.loc[idx_test_positive_SS,'result'] = 1
        self.data_test.loc[idx_test_positive_SS,'observed_status'] = 4
        self.data_test.loc[idx_test_positive_SS,'time'] = t         
        self.data_test.loc[idx_test_positive_SS,'type'] = 0
                
        idx_test_positive_SM = self.data_test.ID.isin(individuals_SM)        
        self.data_test.loc[idx_test_positive_SM,'pos_tim_inf'] = [max(1 + self.data_test["time"][i], t - self.day_since_infected_sympthoms_SM) for i in individuals_SM]
        self.data_test.loc[idx_test_positive_SM,'result'] = 1
        self.data_test.loc[idx_test_positive_SM,'observed_status'] = 5
        self.data_test.loc[idx_test_positive_SM,'time'] = t
        self.data_test.loc[idx_test_positive_SM,'type'] = 0
        
        return True
    
    def test_observations_house_upgrade(self, TP_house, TN_house, status, t):
                 
        idx_test_positive = self.data_test.ID.isin(TP_house)       
        houses_infected = self.data_test.loc[idx_test_positive,'house_no']
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = [max(self.data_test[self.data_test['house_no'] == i]['pos_tim_inf']) for i in np.sort(houses_infected)]
        self.data_test.loc[idx_test_positive,'result'] = 1
        self.data_test.loc[idx_test_positive,'observed_status'] = status[TP_house]
        self.data_test.loc[idx_test_positive,'time'] = t
        self.data_test.loc[idx_test_positive,'type'] = 2 
         
        idx_test_negative = self.data_test.ID.isin(TN_house)        
        self.data_test.loc[idx_test_negative,'result'] = 0
        self.data_test.loc[idx_test_negative,'time'] = t
        
        return True
    
    def test_observations_risk_upgrade(self, TP_risk, TN_risk, status, t ):
             
        idx_test_positive = self.data_test.ID.isin(TP_risk)         
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = [min(self.first_interactions[self.first_interactions['ID_2'] == i]['time'], default = max(0, t - self.gamma )) for i in np.sort(TP_risk)]
        self.data_test.loc[idx_test_positive,'result'] = 1
        self.data_test.loc[idx_test_positive,'observed_status'] = status[TP_risk]
        self.data_test.loc[idx_test_positive,'time'] = t
        self.data_test.loc[idx_test_positive,'type'] = 1 
           
        idx_test_negative = self.data_test.ID.isin(TN_risk)        
        self.data_test.loc[idx_test_negative,'result'] = 0
        self.data_test.loc[idx_test_negative,'time'] = t
       
        return True
    
    def removals_upgrade(self,  status, t ):
        
        all_TP = np.array(self.data_test[self.data_test["result"] == 1]["ID"])
        individuals_removed_at_t = list(np.setdiff1d(all_TP[status[all_TP] >= 8], self.individuals_removed))
        idx = self.data_test.ID.isin(individuals_removed_at_t)
        self.data_test.loc[idx, "time_removed"] = t    
        self.individuals_removed += individuals_removed_at_t
        
        individuals_sucept = np.array(self.data_test[self.data_test["time_infected"] == -1]["ID"])
        individuals_infected_at_t = individuals_sucept[status[individuals_sucept] > 0]
        idx = self.data_test.ID.isin(individuals_infected_at_t)
        self.data_test.loc[idx, "time_infected"] = t    
        
        return len(self.individuals_removed)
    
    def real_time_infection_upgrade(self, TP):
        
        idx_test_positive = self.data_test.ID.isin(TP) 
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = self.data_test.loc[idx_test_positive,'time_infected']
        
        return True
    
 
    def constant_time_infection_upgrade(self, t, TP_risk):
        
        idx_test_positive = self.data_test.ID.isin(TP_risk) 
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = max(0, t - self.gamma)
        
        return True
    
    def rank(self, t, daily_contacts):
        '''
        to compute the 1°CT risk of individuals 
        return: list of ranked individuals
        '''
        
        #To save recent interactions
        self.recent_contacts = [(col[0], col[1], col[2], col[3]) for col in self.recent_contacts if col[2] != t - self.gamma - 1]     
        self.recent_contacts.extend(daily_contacts)
        
        if t >= self.t_start:            
            
            self.first_interactions = f_first_interactions(self.recent_contacts, self.data_test, t, self.prob_function, self.Sa)             
            df_Risk = f_first_contacts_risk(self.first_interactions, self.data_test)
           
            df_Risk['noisse'] = np.random.random(1) - self.noisse
            df_Risk.loc[df_Risk['noisse'] <0, 'noisse'] = 2 + df_Risk.loc[df_Risk['noisse']<0, 'noisse']
            data_to_test = df_Risk[df_Risk['Total_Risk'] >= 0]                
            idx = np.lexsort((data_to_test['noisse'], data_to_test['Total_Risk']))
            individuals_rank = list(data_to_test.iloc[idx]['ID'])
            individuals_rank.reverse()
                        
        else:
            
            individuals_rank = []  
            df_Risk = []
            
        return  individuals_rank
                
    def save_results(self, timeseries, t_end):        
                          
        name_file_res = 'timeseries_1_CT_p_'+str(self.prob_function)+'_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_ti_'+str(self.time_infection)+'_qri_'+str(self.quarantined_random_interactions)
        name_data_test = 'data_test_1_CT_p_'+str(self.prob_function)+'_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_ti_'+str(self.time_infection)+'_qri_'+str(self.quarantined_random_interactions)
        
        df_timeseries = pd.DataFrame(timeseries)
        
        if len(df_timeseries) < t_end:
    
            df_timeseries = pd.merge(pd.DataFrame({'time': range(t_end)}), df_timeseries, how="outer")
        
        df_timeseries.to_csv(self.output_dir + name_file_res+".csv")        
        self.data_test.to_csv(self.output_dir + name_data_test+".csv")        
        
        return df_timeseries
            
    def get_info(self):
        
        return self.data_test, self.first_interactions
    
    def name(self):
        
        return '1°CT'