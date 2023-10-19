#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:53:06 2023

@author: gabriela
"""

import numpy as np
import pandas as pd
from.utils import f_compute_transmission_probability, f_time_to_detection, f_compute_Sa
from.First_Degree_Contact_Tracing_Risk import f_first_contacts_risk, f_first_interactions

def f_risk_first_contacts_by_time(df_FI):
    """
    This function compute:
      - the probability of being infected for a 1°contact at each time 
      
    inputs:
        df_FI : first degree interactions dataframe 
        
    outputs:
        df_FI : data frame for first degree contacts, with the probability of being infected at each time 
    """
    
    df_FI = df_FI.sort_values(['ID_2', 'time'])
    df = df_FI.groupby(['ID_2', 'time']).agg({'CProb': 'prod'})
    
    df_FI =  df_FI.drop_duplicates(subset =["time", 'ID_2'],
                     keep = 'last', inplace = False)  
   
    df_FI['Prob'] = np.array(1 - df['CProb'])
    df_FI['CProb'] = np.array(df['CProb'])
    
    return df_FI


def f_CRisk_1(df_FI):
    """
    This function compute:
      - the cumulative product of 1-probability of being infected for each first case in time 
      
    inputs:
        df_FI : first degree interactions dataframe 
        
    outputs:
        df_Comp_Risk_1 : data frame with columns
                         ID_2: column of 1°contacts 
                         time: time
                         Prob: probability of being infected at each time
                         Comp_Risk_1: the cumulative product of 1-probability of being infected for each 1°contact at each time
                         D_Comp_Risk_1: Comp_Risk_1 outdated one day and putting 1 in the first position
    """
    
    df_FC = f_risk_first_contacts_by_time(df_FI)   
    df_FC['Comp_Risk_1'] = df_FC.groupby(['ID_2']).agg({'CProb': 'cumprod'})
   
    df = df_FC.groupby(['ID_2'])
    df = df['Comp_Risk_1'].apply(lambda x: np. concatenate( ([1], x[:len(x)-1]) )).to_frame().reset_index()
   
    v = df['Comp_Risk_1']
    df_FC['D_Comp_Risk_1'] = np.concatenate(v)

    df_Comp_Risk_1 = df_FC[['ID_2', 'time', 'Prob', 'Comp_Risk_1', 'D_Comp_Risk_1']]
    
    return df_Comp_Risk_1



def f_second_interactions( df_FI_W , Data, df_test, t, zeta):
         
    """
    This function compute the 2°interactions,
         
    inputs:
        df_FI_W  : 1°interactions
        Data : all interactions from the Open ABM Model
        df_test : Data frame with test information
        t : day of intervention
        zeta: time-frame for the time of infection for the 1°contact
        
    outputs:
        df_SI_W : The second degree intercations 
        
    """
    Data = pd.DataFrame(Data, columns = ['ID', 'ID_2', 'time', 'type'])   
    
    # all tested positives 
    all_TP = np.array(df_test[(df_test['result'] == 1)]['ID'])     
      
    if len(df_FI_W)>0:
        
        #to filter the second degree interactions
        time_active = list(range(t - zeta , t + 1))        
        fc = df_FI_W ['ID_2'].unique() #first degree contacts    
        df_SI = Data[Data.time.isin(time_active) & Data.ID.isin(fc) & ~Data.ID_2.isin(all_TP) ]     
        df_SI = df_SI.sort_values(['ID']) # sort the data frame by the first degree contacts
            
        if len(df_SI) > 0 :
            
            # to filter the interactions after the test negative day for second degree contacts            
            df_SI = df_SI.merge(df_test[['ID', 'time']].rename(columns = { 'ID': 'ID_2', 'time':'time_test_neg_2'}))
            df_SI = df_SI[df_SI['time_test_neg_2'] < df_SI ['time']] 
               
            if len(df_SI) > 0:
                
                #to consider asymtomatics status for each first degree contact
                df_SI['observed_status_1'] = np.zeros(len(df_SI))  
                df2 = df_test[['ID', 'age_group']].rename(columns = {'ID': 'ID_2', 'age_group':'age_group_2'})
                df_SI = df_SI.merge(df2)
                                
    else:
        
        df_SI = pd.DataFrame()          
    
    return df_SI

def f_second_risk_data(df_FI, df_SI_W, prob_function, Sa):
     
    """
    This function compute the 2°CT risk formula for each interaction between a 1° and a 2° contact,         
    inputs:
        df_FI_W : 1°interactions
        df_SI_W : 2°interactions         
    outputs:
        data : the Column 'Comp_Risk_2_1' represent the probability that each 1° contact  non infects each 2° contact in all posibles interactions times
             
    """ 
    
    if prob_function == 0:
        prob_function = -1
    
    SI = df_SI_W[['ID', 'ID_2', 'time', 'type', 'observed_status_1', 'age_group_2']]
   
    FI = f_CRisk_1(df_FI)
    FI.columns = ['ID', 'time_1', 'Prob1', 'CRisk_1', 'DCR_1']
    
    data = SI.merge(FI)
    data = data.sort_values(['ID_2','ID', 'time_1', 'time'])
    data = data[data['time_1'] < data['time']]
    
    data['pos_tim_inf_1'] = data['time_1']
        
    data['Prob2'] =  f_compute_transmission_probability(data, prob_function, Sa)   
    data['CProb2'] = 1 - data['Prob2']
     
    data = data.sort_values(['ID_2','ID', 'time_1', 'time'])
    
    df_r2 = data.groupby(['ID_2','ID', 'time_1'], as_index=False).agg({'CProb2': 'prod', 'Prob1': 'min', 'DCR_1': 'min'})
    
    df_r2['Mult'] = df_r2['Prob1'] * df_r2['DCR_1']* df_r2['CProb2'] 
    
    df_risk = df_r2.groupby(['ID_2','ID'], as_index=False).agg({'Mult': 'sum', 'time_1': 'max'})
    
    df_risk = df_risk.merge(FI[['ID', 'time_1', 'CRisk_1']])
    df_risk['Comp_Risk_2_1'] = df_risk['Mult'] + df_risk['CRisk_1'] 
    
    return df_risk

def f_second_contacts_risk(df_FI, df_SI_W, n_total, prob_function, Sa):
    
    """
    This function compute the risk for 2°contacts
    
    """    
    data = f_second_risk_data(df_FI, df_SI_W, prob_function, Sa)
   
    df = data.groupby(["ID_2"], as_index=False).agg({'Comp_Risk_2_1': 'prod'})    
    df.columns = ['ID', 'Risk_2']  
    df['Risk_2'] = 1 - df['Risk_2']
    
    # Put zero in the nan values
    df_Risk_2 = pd.merge(pd.DataFrame({'ID': range(n_total)})  , df, how="outer")
    df_Risk_2['Risk_2'] = df_Risk_2['Risk_2'].fillna(0)
    
    return df_Risk_2


def f_global_risk(df_FI, df_SI, data_test, prob_function, Sa, t, gamma):
    
    """
    This function compute the risk for 1° and 2° contacts
    
    """
    
    n_total = len(data_test['ID'])
    
    if len(df_FI) > 0:
        df_Risk_1 = f_first_contacts_risk(df_FI[df_FI['time'] >= max(0, t - gamma)], data_test)
    else:
        df_Risk_1 = pd.DataFrame({'ID': range(n_total), 'Total_Risk': 0})
        
    if len(df_SI) > 0:    
        df_Risk_2 = f_second_contacts_risk( df_FI, df_SI, n_total, prob_function, Sa)   
    else:
        df_Risk_2 = pd.DataFrame({'ID': range(n_total), 'Risk_2': 0})   
    
    #Second degree contacts risk        
    df_Risk = df_Risk_2        
    df_Risk['Risk_1'] =  df_Risk_1['Total_Risk']  

    #Second degree contact tracing risk      
    df_Risk['Total_Risk'] = 1 - (1 - df_Risk['Risk_1'])*(1 - df_Risk['Risk_2'])    
        
    all_TP = data_test[data_test["result"] == 1]["ID"]    
    df_Risk.loc[df_Risk.ID.isin(all_TP), 'Total_Risk'] = -1  
    
    return df_Risk


class Second_Degree_Contact_Tracing_Risk:
        
    def __init__(self, gamma, zeta, prob_function):
        self.description = "class for second degree contact tracing method."   
        self.gamma = gamma 
        self.zeta = zeta        
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
        self.recent_contacts = []
        self.first_interactions = []
        self.recent_first_interactions = []
        self.second_interactions = []
        self.individuals_removed = list()
        self.rng = np.random.RandomState(1)
        
        if self.prob_function == 'lambda':
           self.Sa = f_compute_Sa(age_group, house_no)
           
        else:
            self.Sa = []
        
        data_test = pd.DataFrame()
        data_test["ID"] = np.arange(n_total)
        data_test["result"] = np.zeros(n_total) - 1
        data_test["time"] = np.zeros(n_total) - 1
        data_test["pos_tim_inf"] = np.zeros(n_total) - 1
        data_test["observed_status"] = np.zeros(n_total) 
        data_test["age_group"] = age_group
        data_test["house_no"] = house_no
        data_test["time_removed"] = np.zeros(n_total) + np.inf
        data_test["time_infected"] = np.zeros(n_total) - 1
        data_test["type"] = np.zeros(n_total) - 1 #0: detection by symptomps, 1: by risk, 2: by home, -1: not detected
        self.data_test = data_test
        
        np.random.seed(1)
        self.noisse = np.random.random(n_total)
       
        return True
    
    def test_observations_symptomps_upgrade(self, individuals_SS, individuals_SM, t):
        
        idx_test_positive_SS = self.data_test.ID.isin(individuals_SS)        
        self.data_test.loc[idx_test_positive_SS,'pos_tim_inf'] = [max(1 + self.data_test["time"][i], t - self.day_since_infected_sympthoms) for i in individuals_SS]
        self.data_test.loc[idx_test_positive_SS,'result'] = 1
        self.data_test.loc[idx_test_positive_SS,'observed_status'] = 4
        self.data_test.loc[idx_test_positive_SS,'time'] = t
        self.data_test.loc[idx_test_positive_SS,'type'] = 0
        
        idx_test_positive_SM = self.data_test.ID.isin(individuals_SM)        
        self.data_test.loc[idx_test_positive_SM,'pos_tim_inf'] = [max(1 + self.data_test["time"][i], t - self.day_since_infected_sympthoms) for i in individuals_SM]
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
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = [min(self.recent_first_interactions[self.recent_first_interactions['ID_2'] == i]['time'], default= np.inf) if len(self.recent_first_interactions[self.recent_first_interactions['ID_2'] == i])>0 else min(self.second_interactions[self.second_interactions['ID_2'] == i]['time'], default = max(0, t - self.gamma)) for i in np.sort(TP_risk)]
        self.data_test.loc[idx_test_positive,'result'] = 1
        self.data_test.loc[idx_test_positive,'observed_status'] = status[TP_risk]
        self.data_test.loc[idx_test_positive,'time'] = t
        self.data_test.loc[idx_test_positive,'type'] = 1 
            
        idx_test_negative = self.data_test.ID.isin(TN_risk)        
        self.data_test.loc[idx_test_negative,'result'] = 0
        self.data_test.loc[idx_test_negative,'time'] = t
       
        return True
    
    def removals_upgrade(self, status, t ):
               
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
        self.data_test.loc[idx_test_positive,'pos_tim_inf'] = max( 0, t - self.gamma)
        
        return True
    
    def rank(self, t, daily_contacts):
        '''
        to compute the 2CT risk of individuals 
        return: list of ranked individuals
        '''
        
        #To save recent interactions
        self.recent_contacts = [(col[0], col[1], col[2], col[3]) for col in self.recent_contacts if col[2] != t - self.zeta -1]     
        self.recent_contacts.extend(daily_contacts)
        
        if t >= self.t_start:            
            
            #First degree interactions    
            self.first_interactions  = f_first_interactions(self.recent_contacts, self.data_test, t , self.prob_function, self.Sa)             
            
            if len(self.first_interactions)>0:
                self.recent_first_interactions = self.first_interactions[self.first_interactions['time'] >= t - self.gamma]
            else :
                self.recent_first_interactions = self.first_interactions
                
            #Second degree interactions         
            self.second_interactions = f_second_interactions(self.first_interactions, self.recent_contacts, self.data_test, t, self.gamma)  
            
            #to compute the risk            
            df_Risk = f_global_risk(self.first_interactions, self.second_interactions, self.data_test, self.prob_function, self.Sa, t, self.gamma)
            
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
                          
        name_file_res = 'timeseries_2_CT_p_'+str(self.prob_function)+'_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_zeta_'+str(self.zeta)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_ti_'+str(self.time_infection)+'_qri_'+str(self.quarantined_random_interactions)                  
        name_data_test = 'data_test_2_CT_p_'+str(self.prob_function)+'_n_total_'+str(self.n_total)+'_N0_'+str(self.initial_inf)+'_t_start_'+str(self.t_start)+'_eta_'+str(self.test_availables)+'_gamma_'+str(self.gamma)+'_zeta_'+str(self.zeta)+'_qh_'+str(self.quarantine_household)+'_th_'+str(self.test_household)+'_p_SS_'+str(self.p_SS)+'_p_SM_'+str(self.p_SM)+'_seed_'+str(self.seed)+'_seed_open_ABM_'+str(self.seed_open_ABM)+'_ti_'+str(self.time_infection)+'_qri_'+str(self.quarantined_random_interactions)                  
                       
        df_timeseries = f_time_to_detection(timeseries, self.data_test)    
        df_timeseries['det_prop'] = (df_timeseries['detected_R'] + df_timeseries['detected_H'])/(df_timeseries['active'] - df_timeseries['detected_active'] + df_timeseries['detected_R']+ df_timeseries['detected_H'])
    
        if len(df_timeseries) < t_end:
    
            df_timeseries = pd.merge(pd.DataFrame({'time': range(t_end)}), df_timeseries, how="outer")
        
        df_timeseries.to_csv(self.output_dir + name_file_res+".csv")        
        self.data_test.to_csv(self.output_dir + name_data_test+".csv")        
        
        return df_timeseries
      
    def name(self):
        
        return '2°CT'     
    
    def get_info(self):
        
        return self.first_interactions, self.second_interactions, self.data_test        
  