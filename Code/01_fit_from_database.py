from Thermal_Models import estimate_parameters, \
                           schoolfield_two_factor, \
                           schoolfield_two_factor_tpk, \
                           schoolfield_full, \
                           schoolfield_full_tpk, \
                           read_database, \
                           fit_models, \
                           split_datasets, \
                           rank_and_flatten, \
                           compile_models, \
                           estimate_uncertainty, \
                           bootstrap_model
import numpy as np
import pandas as pd
from datetime import datetime
import warnings

warnings.simplefilter(action = "ignore", category = FutureWarning) #Probably ought to deal with this, but for the moment I'll sweep it under the rug. 
warnings.simplefilter(action = "ignore", category = RuntimeWarning )

starttime = datetime.now()

data_path = '../Data/database.csv'
data = read_database(data_path) 
#data = data[data.StandardisedTraitName == 'Specific Growth Rate']
Datasets = split_datasets(data) #return a dictionary of growth curves split by ID

#names of the models I want to fit as strings
model_names = ['schoolfield_two_factor_tpk']
all_models = [] #place to put fitted models

#Column names of parameters I want to use as explanatory variables when I analyse the data                  
aux_parameters = ['FinalID', 'OriginalID', 'Citation', 'Latitude', 'Longitude', 'ConKingdom', 'ConPhylum', 'ConClass',
                  'ConOrder', 'ConFamily', 'ConGenus', 'ConSpecies', 'OptimalConditions', 'RespiratoryPathway', 'ConLabGrowthTemp'] 
                  
for i in Datasets.keys(): #Iterate through the dictionary by key
    dataset = Datasets[i] #get the growth curve
    est_params = estimate_parameters(dataset, aux_parameters) #Estimate starting parameters using regression
    models = fit_models(model_names, est_params, use_randomisation=True, tag = i) #Fit 3 models in one line, returns a list containing them
    
    if models: #Check something has actually fitted
        best_model = max(models) #Find the best model
        models = [estimate_uncertainty(best_model, est_params, bootstrap_allowed = False)]
        
        if best_model.final_E > 0.001: #Checks that the model actually has a growth response
            all_models.append(models) #could also be all_models.append(best_model) if you aren't interested in the alternate fits
        print(best_model)
        #best_model.plot('../Results/fits', convert_kelvin = True, show_estimates=False) # want to see the fits in celsius, though currently incompatible with plotting uncertainty - TS
         
all_models = rank_and_flatten(all_models)
output_data = compile_models(all_models, show_estimates=True, aux_cols = aux_parameters) #Save the summary in two places

output_data.to_csv('../Results/summary.csv')
output_data.to_csv('../Data/summaries/summary.csv')

print('Completed in: ', datetime.now() - starttime)
