#!/usr/bin/env python3
"""
This program provides a framework for and implementation of least squares fitting of various 
thermal response models based on experimental data

Written in Python 3.5 Anaconda Distribution
"""

from Thermal_Models import estimate_parameters, \
                           read_database, \
                           fit_models, \
                           bootstrap_model, \
                           split_datasets, \
                           rank_and_flatten, \
                           compile_models, \
                           find_linear_arrhenius, \
                           estimate_uncertainty
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

starttime = datetime.now()

summary_path = '../Data/summaries/group_means_non_aggregated.csv' #We will merge results into this
data_path = '../Data/summaries/non_aggregate_data_fix.csv'
data = read_database(data_path) 
summary = read_database(summary_path) 
data = data[(np.isfinite(data.Est_Tpk)) & (np.isfinite(data.Est_Response))]

aux_parameters = [] #no supplementary information needed
analysis_levels = ['ConKingdom'] #levels at which models should be fitted 
model_names = ['boltzmann_arrhenius_one_weight_non_linear'] #models to fit

for level in analysis_levels:
    if level == '':
        hist_axes = True
    else:
        hist_axes = False
    
    #see readme for more details on this dictionary
    param_est_flags = {'Trait_Col_Name': 'Trait',
                       'X_vals_col_name': 'Est_Tpk',
                       'Y_vals_col_name': 'Est_Response',
                       'x_val_conversion': 1,
                       'log_data': False,
                       'is_celcius': False,
                       'species_data': False,
                       'full_strain_col_name': level,
                       'genus_level_col_name': level,
                       'species_level_col_name': level,
                       'uncertainty_x_col_name': 'Est_Tpk_std',
                       'uncertainty_y_col_name': 'Est_Response_std'
                       }
                   
    all_models = []
    
    Datasets = split_datasets(data, sep = level, _sort = ['Est_Tpk'])      
    
    count = 1
    
    for i in Datasets.keys():
        dataset = Datasets[i]
        tmax_check = len(dataset['Est_Tpk'].dropna())
        
        if tmax_check >= 5: #Must have more datapoints than number of variables
            est_params = estimate_parameters(dataset, aux_parameters, flags = param_est_flags)
            #est_params = find_linear_arrhenius(est_params, show=False)
             
            #print(est_params)
            
            models = fit_models(model_names, est_params, tag = (i + 1))
            models = [bootstrap_model(model, est_params, N=1000) for model in models] #boostrap each model
            
            all_models.append(models)
            
            plot_path = '../Results/without_aggregation/{}/standard'.format(level)
            plot_path_ahhr = '../Results/without_aggregation/{}/arrh'.format(level)
            
            for i, model in enumerate(models):
                model.plot(plot_path, plot_residuals=False, hist_axes = hist_axes, fit_stats = True, convert_kelvin = True, plot_value_uncertainty = True, manual_axis_limits=True)
                model.plot(plot_path_ahhr, plot_residuals=False, hist_axes = hist_axes, fit_stats = True, convert_kelvin = False, scale_type = 'arrhenius', manual_axis_limits=False)

                
            count += 1

    all_models = rank_and_flatten(all_models)
    output_data = compile_models(all_models, aux_cols = aux_parameters, sortby=['Species', 'Model_name'])
    
    summary_path_extra = '../Results/Maxima_fits/{}_summary_BA_OWNL.csv'.format(level)
    output_data.to_csv(summary_path_extra)
    
    #select appropriate columns
    merge_data = output_data[['Species', 'Model_name', 'E', 'E_min', 'E_max', 'B0']]
    merge_data.rename(columns={'Species': level,
                               'Model_name': 'Model_name',
                               'E': level + '_EG',
                               'E_min': level + '_EG_min',
                               'E_max': level + '_EG_max',
                               'B0': level + '_B0'}, inplace=True)
 
    merge_data.replace('Deinococcus-thermus', 'Deinococcus-Thermus', inplace=True) #lazy fix
    
    #merge into summary - use this for multiple models
    #~ if level == 'ConKingdom':
        #~ summary = pd.merge(summary, merge_data, how='outer', on=[level])
    #~ else:
        #~ pd.merge(summary, merge_data, how='outer', on=[level, 'Model_name'])
    
    #merge into summary - use this for one model
    summary = pd.merge(summary, merge_data, how='outer', on=level)
    
#split by temppref
Datasets = split_datasets(data, sep = 'ConKingdom', _sort = ['Est_Tpk'])  
TempPref_ach = split_datasets(Datasets[1], sep = 'TempPref', _sort = ['Est_Tpk'])
TempPref_bac = split_datasets(Datasets[0], sep = 'TempPref', _sort = ['Est_Tpk'])    
all_temps = [('{}_archaea'.format(i), TempPref_ach[i]) for i in TempPref_ach.keys()] + [('{}_bacteria'.format(i), TempPref_bac[i]) for i in TempPref_bac.keys()]
count = 0  
temp_models = []

for group in all_temps:
    count += 1
    param_est_flags = {'Trait_Col_Name': 'Trait',
                       'X_vals_col_name': 'Est_Tpk',
                       'Y_vals_col_name': 'Est_Response',
                       'x_val_conversion': 1,
                       'log_data': False,
                       'is_celcius': False,
                       'species_data': False,
                       'full_strain_col_name': 'TempPref',
                       'genus_level_col_name': 'ConKingdom',
                       'species_level_col_name': 'TempPref',
                       'uncertainty_x_col_name': 'Est_Tpk_std',
                       'uncertainty_y_col_name': 'Est_Response_std'
                       }

    data = group[1]
    key = group[0]
    tmax_check = len(data['Est_Tpk'].dropna())
    
    if tmax_check >= 5: #Must have more datapoints than number of variables
        est_params = estimate_parameters(data, flags = param_est_flags, aux_parameters_names=['ConKingdom'])
        #est_params = find_linear_arrhenius(est_params, show=False) 
        #print(est_params)
        
        models = fit_models(model_names, est_params, tag = count)
        models = [bootstrap_model(model, est_params, N=1000) for model in models]
        temp_models.append(models)
        
        plot_path = '../Results/without_aggregation/TempPref/standard'
        plot_path_ahhr = '../Results/without_aggregation/TempPref/arrh'
            
        for i, model in enumerate(models):
            model.plot(plot_path, plot_residuals=False, hist_axes = hist_axes, fit_stats = True, convert_kelvin = True, plot_value_uncertainty = True, manual_axis_limits=True)
            model.plot(plot_path_ahhr, plot_residuals=False, hist_axes = hist_axes, fit_stats = True, convert_kelvin = False, scale_type = 'arrhenius', manual_axis_limits=False)
        
temp_models = rank_and_flatten(temp_models)
output_data = compile_models(temp_models, sortby=['Species', 'Model_name'], aux_cols=['ConKingdom']) 
merge_data = output_data[['Species', 'Model_name', 'ConKingdom', 'E', 'E_min', 'E_max', 'B0']]

merge_data.rename(columns={'Species': 'TempPref',
                           'Model_name': 'Model_name',
                           'ConKingdom': 'ConKingdom', 
                           'E': 'TempPref_EG',
                           'E_min': 'TempPref_EG_min',
                           'E_max': 'TempPref_EG_max',
                           'B0': 'TempPref_B0'}, inplace=True)
                           
summary = pd.merge(summary, merge_data, how='outer', on=['TempPref', 'ConKingdom'])



summary.fillna('NA', inplace=True)    
summary.replace(r'\s+', 'NA', regex=True, inplace=True)
summary.to_csv('../Data/summaries/all_means_non_aggregated.csv',  index = False)

print('Completed in: ', datetime.now() - starttime)
