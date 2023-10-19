#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 17:26:04 2023

@author: gabriela
"""

import pandas as pd
import numpy as np
from scipy.stats import gamma
import scipy.sparse  

def f_compute_Sa(age_group, house_no):
        
    df_indiv = pd.DataFrame({'ID': np.arange(len(age_group)), 'age_group': age_group, 'house_no': house_no})     
    
    df_indiv['new_age_group'] = np.where(df_indiv.age_group.isin([0,1]) , 0, np.where(df_indiv.age_group.isin([2,3,4,5,6]) , 1, 2))
    df_house = df_indiv[['ID', 'house_no']].groupby('house_no', as_index=False).count()
    df_house['house_interactions'] = df_house['ID'] - 1
    df = df_indiv.merge(df_house[['house_no', 'house_interactions']])
    mean = df[['new_age_group', 'house_interactions']].groupby('new_age_group', as_index=False).mean()
    mean['total_interactions'] = np.array([12, 11, 6])
    mean['total_interactions'] = mean['house_interactions'] + mean['total_interactions']
    I = np.repeat(np.array(mean['total_interactions']), [2, 5 , 2])
    
    S = [0.71, 0.74, 0.79, 0.87, 0.98, 1.11, 1.26, 1.45, 1.66]   
    Sa = np.divide(S, I)        
    
    return Sa

def f_compute_weigth_ABM(time, status, age, network, Sa):
    
    "Open ABM transmission probability"
    
    prob = []    
    f =[]    
    time = np.array(time)    
    n_rows = len(status)        
    # relative infectiousness of the source based on the disease severity (asymptomatic, mild, moderate/severe)
    # 0 (suceptible), 3 (asymptomatic), 1-4 (pre-severe,severe), 2-5 (pre-mild,mild) 
    
    k = 1    
    #Ad = [0.25 * k, k, 0.48* k, 0.25* k, k, 0.48 * k, 0, 0, 0, 0, 0]
    Ad = [0.29 * k, k, 0.48* k, 0.29* k, k, 0.48 * k, 0, 0, 0, 0, 0]
    
    matrix_status = np.zeros((n_rows, len(Ad)))    
    status = np.array(status)    
    data = np.ones(n_rows)    
    row = np.arange(n_rows)    
    status_sparse_matrix = scipy.sparse.coo_matrix((data, (row, status)), 
                                                           shape=matrix_status.shape, 
                                                           dtype=np.int8)  
    vector_status = status_sparse_matrix * Ad
    
    #relative susceptibility of the recipient base on age    
    matrix_age = np.zeros((n_rows, len(Sa)))    
    age = np.array(age)    
    data = np.ones(n_rows)    
    row = np.arange(n_rows)        
    age_sparse_matrix = scipy.sparse.coo_matrix((data, (row, age)), 
                                                        shape=matrix_age.shape, 
                                                        dtype=np.int8)        
    vector_age = age_sparse_matrix * Sa
    
    #scale factor for the network on which the interaction ocurred(n: type of network: random, household and differents types of  works)    
    Bn = [2, 1, 1]    
    matrix_network = np.zeros((n_rows, len(Bn)))    
    network = np.array(network)    
    data = np.ones(n_rows)    
    row = np.arange(n_rows)        
    network_sparse_matrix = scipy.sparse.coo_matrix((data, (row, network)), 
                                                           shape=matrix_network.shape, 
                                                           dtype=np.int8)        
    vector_network = network_sparse_matrix * Bn     
    
    #scales the overal infectious rate
    R = 5.75              
    
    #mean and width parameters from gamma density function
    mn = 6    
    sd = 2.5      
    scale = (sd**2)/mn
    shape = mn/scale    
    #f = gamma.cdf(time , shape, scale) - gamma.cdf(time - 1, shape, scale) 
    f = gamma.cdf(time + 1 , shape, scale) - gamma.cdf(time, shape, scale)                   
    lambd = R*(vector_age)*(vector_status)*(vector_network)*f
    
    prob = 1 - np.exp(-lambd)
    
    return prob

def f_compute_transmission_probability(df, prob_function, Sa):
    
    if prob_function == 'cte':
        
       result = 1/2   
       prob =  pd.Series(result, index=df.index, name='result')
     
    if prob_function == 'lambda':
      
       result = f_compute_weigth_ABM(df['time'] - df['pos_tim_inf_1'], df ['observed_status_1'], df['age_group_2'], df['type'], Sa)   
       prob =  pd.Series(result, index=df.index, name='result')   
    
    return prob
    
def f_time_to_detection(data_timeseries, data_test):
    
    timeseries_sim = pd.DataFrame(data_timeseries)
    
    t_end = len(timeseries_sim)    
    data_test['time_to_detection'] =  data_test['time']  - data_test['time_infected']
    
    #individuals detected
    data = data_test[data_test['result'] == 1][['time', 'time_to_detection']]       
    df = data.groupby(['time'], as_index=False).agg('sum')
    df.columns = ['time', 'sum_time']        
    
    # Put nan in others values
    df_new = pd.merge(pd.DataFrame({'time': range(t_end)}), df, how="outer")
    timeseries_sim['sum_time_det'] = df_new['sum_time']
        
    #individuals detected by symptoms
    data_S = data_test[data_test['type'] == 0][['time', 'time_to_detection']]       
    #df_S = data_S.groupby(['time'], as_index=False).agg({'time_to_detection': 'mean'})
    df_S = data_S.groupby(['time'], as_index=False).agg('sum')
    df_S.columns = ['time', 'sum_time']
       
    # Put nan in others values
    df_Symp = pd.merge(pd.DataFrame({'time': range(t_end)})  , df_S, how="outer")
    #df_Symp.fillna(0)
    timeseries_sim['sum_time_det_S'] = df_Symp['sum_time']
    
    #individuals detected by risk
    data_R = data_test[data_test['type'] == 1][['time', 'time_to_detection']]       
    df_R = data_R.groupby(['time'], as_index=False).agg('sum')
    df_R.columns = ['time', 'sum_time']
       
    # Put nan in others values
    df_Risk = pd.merge(pd.DataFrame({'time': range(t_end)})  , df_R, how="outer")
    #df_Symp.fillna(0)
    timeseries_sim['sum_time_det_R'] = df_Risk['sum_time']
    
    return timeseries_sim