import numpy as np
import pandas as pd
import itertools
import warnings

from progress.bar import Bar
from copy import deepcopy
from scipy import stats

from .models import boltzmann_arrhenius, lm, schoolfield_two_factor, schoolfield_two_factor_tpk, schoolfield_full, schoolfield_full_tpk, boltzmann_arrhenius_two_weights_linear, boltzmann_arrhenius_two_weights_non_linear, boltzmann_arrhenius_two_weights_normalise_only, boltzmann_arrhenius_one_weight_linear, boltzmann_arrhenius_one_weight_non_linear

import matplotlib.pyplot as plt
import seaborn as sns


plotname = 1

#Dictionary of model to fit, min number of variables
all_models = {
            'lm' : (lm, 2),
            'boltzmann_arrhenius' : (boltzmann_arrhenius, 2),
            'schoolfield_two_factor_tpk': (schoolfield_two_factor_tpk, 4),
            'schoolfield_two_factor': (schoolfield_two_factor, 4), 
            'schoolfield_full': (schoolfield_full, 6), 
            'schoolfield_full_tpk': (schoolfield_full_tpk, 6),
            'boltzmann_arrhenius_two_weights': (boltzmann_arrhenius_two_weights_non_linear, 2), #redundancy entry for backwards compatibility
            'boltzmann_arrhenius_two_weights_non_linear': (boltzmann_arrhenius_two_weights_non_linear, 2),
            'boltzmann_arrhenius_two_weights_linear': (boltzmann_arrhenius_two_weights_linear, 2), 
            'boltzmann_arrhenius_two_weights_normalise_only': (boltzmann_arrhenius_two_weights_normalise_only, 2),
            'boltzmann_arrhenius_one_weight_linear': (boltzmann_arrhenius_one_weight_linear, 2),
            'boltzmann_arrhenius_one_weight_non_linear': (boltzmann_arrhenius_one_weight_non_linear, 2)
             }

def read_database(path):
    "Read the file in latin 1, convenience function"
    return pd.read_csv(path, encoding = "ISO-8859-1", low_memory = False)

def fit_models(models, estimates, tag=None, print_each=False, use_randomisation=False, randomisation_tries=50, randomisation_scale=1, prog_bar = True, fit_vars=all_models):
    """
    wrapper for the model classes, performs basic checks and returns a list of fitted models

    models = list/tuple of model names as strings to fit to each dataset
    estimates = list/tuple of parameter estimate objects for each dataset
    tag = A tag to use as a name for any plots produced
    print_each = bool flag to determine if fit details should be printed to console
    
    """
    assert type(models) in (list, tuple), "Models must be passed as a list of strings!"
    
    models = [i.lower() for i in models] #convert model names to lower case for robustness
    fitted_models = [] #Use to store model objects
    not_fitted = [] #List to detail any models which fail to fit
    n_vars = len(estimates.temps) #Number of data points in the dataset
    
    if tag == None: #Make sure plot has a name
        tag = "{}_".format(estimates.trait)

    if use_randomisation:
        random_params = randomise_estimates(estimates, N=randomisation_tries, scale_levels=randomisation_scale)
    
    for i in models:
        model_class, min_vars = fit_vars[i]
        if n_vars >= min_vars: #check there are enough data points to fit the model
            model = model_class() #initiate the model, don't do anything
            model.fit_from_parameters(estimates, tag) #perform the fit

            if use_randomisation and ((model.R2 < 0.8) or not (model.model.covar is not None and np.all(np.isfinite(model.model.covar)))):
                print('\n\tTrying random parameters!')

                if prog_bar:
                    progbar = Bar('\tTrying to improve model!', max=randomisation_tries)

                best_model = model

                j = 0 #incrementor
                covar_mat = best_model.model.covar

                while j < randomisation_tries and ((best_model.R2 < 0.8) or not (covar_mat is not None and np.all(np.isfinite(covar_mat)))):
                    random_parameter_set = random_params[j] #get a random paramter set

                    random_model = model_class() #initiate the model again
                    random_model.fit_from_parameters(random_parameter_set, tag) #fit the model again

                    best_model = max([random_model, best_model]) #check which model is better

                    covar_mat = best_model.model.covar
                    j += 1

                    if prog_bar:
                        progbar.next()

                if prog_bar:
                    progbar.finish()

                model = best_model #find the best fitted model using AIC

            if print_each:
                print(model)

            fitted_models.append(model)
        else: 
            not_fitted.append(i)
            
    if not_fitted: #if list isn't empty print details of failed fits
        not_fitted_str = ', '.join(not_fitted) #list to comma seperated string
        print('\tNot enough data points to fit: {}\n'.format(not_fitted_str))
        
    return fitted_models

def set_model_percentiles(model, params_dict):
    "Set limits on each bootstrapped parameter within the model"
    no_boot = set(['NA']) #If the parameter was not boostrapped the list will be entirely NA

    for parameter in params_dict.keys(): #iterate through all parameters passed to function
        param_vals = params_dict[parameter]
        param_max_name = parameter + '_max' #Names for min and max values to be applied to model
        param_min_name = parameter + '_min'

        if set(param_vals) != no_boot: #check the list of values isn't empty
            param_vals = [i for i in param_vals if i not in ('NA', np.inf, -np.inf)]
            upper_lim = np.nanpercentile(param_vals, 97.5) #use nanpercentile for robustness
            lower_lim = np.nanpercentile(param_vals, 2.5)

            #remove this
            #~ if parameter == 'final_E':
                #~ global plotname
                #~ plt.figure()
                #~ plt.hist(param_vals, 30, facecolor='green', alpha=0.75)
                #~ plt.xlabel('E')
                #~ plt.ylabel('Frequency')
                #~ plt.title('Bootstrapped EG value distribution')
                #~ plt.savefig('../Results/EG_Histograms/{}.png'.format(str(model.species_name)))
                #~ plt.close()
                #~ plotname+=1

            setattr(model, param_max_name, upper_lim) #set upper confidence limit
            setattr(model, param_min_name, lower_lim) #set lower confidence limit

    return model 

def perform_bootstrap(model, parameters, N, suppress_progress=False, fit_vars=all_models):
    "Perform a bootstrap on the model, to see the output you need to set bootstrap to true in output_csv"
    
    hashmap = {} #Use this to store completed models, avoid rerunning fits with small datasets to save time
   
    #names of model parameters to bootstrap
    params_to_bootstrap = ['final_E', 'final_B0', 'tpk_est', 'max_response_est', 'final_E_D',
                           'final_E_D_L', 'final_T_H', 'final_T_H_L', 'slope', 'intercept']
    bootstrapped_params = {parameter: list() for parameter in params_to_bootstrap} #create a dict of parameters and empty lists

    #Bootstrapping is slow so use a progress bar to show something is happening!
    #progressbar gets buggy if the text is wider than the output window. 
    bar_lab = '\tBootstrap {} {}'.format(str(model.species_name)[:12], str(model.model_name_short))     

    if not suppress_progress:
        progbar = Bar(bar_lab, max=N)
    else:
        print('\tBootstrapping {} {}\n'.format(str(model.species_name), str(model.model_name_short)))
    
    model_type = model.model_name.lower() #don't want to directly reuse the object as that could have implications elsewhere

    #Iterate N times
    for i in range(N):
        new_parameters = deepcopy(parameters) #copy over parameters as they have built in resample
        new_parameters.resample_data()
        
        unique_key = tuple(np.concatenate((np.sort(new_parameters.temps), np.sort(new_parameters.responses)))) #create an immutable unique key for the hashmap
        
        if len(unique_key) < 30 and unique_key in hashmap: #If we already ran a model with these parameters there is no point running it again
            bootstrap_model = hashmap[unique_key] 

        else: #But otherwise re-estimate everything
            new_parameters.estimate_all(boot=True)
            bootstrap_model = fit_vars[model_type][0]() #Fresh model
            bootstrap_model.fit_from_parameters(new_parameters, index='boot') #Fit the model using new parameters
            hashmap[unique_key] = bootstrap_model  #add the fitted model to the hashmap 
        
        #extract each paramter and add it to a list contained in the boostrapped_params dictionary
        for param in params_to_bootstrap:
            result = getattr(bootstrap_model, param, "NA")
            bootstrapped_params[param].append(result)
        
        if not suppress_progress:
            progbar.next()  
    if not suppress_progress:
        progbar.finish()
    
    return bootstrapped_params

def bootstrap_model(model, parameters, N = 1000, suppress_progress=False):
    """Seperate the boostrap process from the paramter setting process for flexibility"""
    bootstrapped_params = perform_bootstrap(model, parameters, N, suppress_progress=suppress_progress)
    model = set_model_percentiles(model, bootstrapped_params)
    model.bootstrapped = True

    return model

def estimate_uncertainty(model, parameters, bootstrap_allowed = False, n = 1000):
    "Boostrap if the covariance matrix reaches singularity"
    tpk_values, max_response_values, E_values = model.estimate_uncertainty(n)

    tpk_uncertainty = np.nanstd(tpk_values)  #take standard deviation
    response_uncertainty = np.nanstd(max_response_values)
    E_uncertainty = np.nanstd(E_values) 

    if bootstrap_allowed and (not np.isfinite(tpk_uncertainty) or not np.isfinite(response_uncertainty)):
        bootstrapped_params = perform_bootstrap(model, parameters, 1000)

        tpk_uncertainty = np.nanstd(bootstrapped_params['tpk_est'])  #take standard deviation
        response_uncertainty = np.nanstd(bootstrapped_params['max_response_est']) 
        E_uncertainty = np.nanstd(bootstrapped_params['final_E']) 

    model.tpk_uncertainty = tpk_uncertainty
    model.response_uncertainty = response_uncertainty
    model.E_uncertainty = E_uncertainty

    return model
            
def split_datasets(data, sep = 'OriginalID', _sort = ['ConTemp']):
    "Create a set of temperature response curve datasets from csv"
    data['FinalID'] = pd.factorize(data[sep])[0]
    ids = pd.unique(data['FinalID']).tolist() #Get unique identifiers
    
    #create a dictionary of datasets for easy access later
    Datasets = {}
    for id in ids:
        curve_data = data.loc[data['FinalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(_sort).reset_index() #sort so rows are in temperature order, reset index to 0  
        Datasets[id] = curve_data
    return Datasets    

def rank_and_flatten(model_list):
    "Function to rank a nested lit of models and flatten the list"
    all_models = []
    model_list = [i for i in model_list if i != []]
    if any(isinstance(i, list) for i in model_list): #check if list is nested
        for candidates in model_list:
            candidates.sort() #Uses built in lt method, best model will be first
            min_aic = candidates[-1].AIC
            min_bic = candidates[-1].BIC
            best_model_type = candidates[-1].model_name
            for rank, model in enumerate(candidates[::-1]):
                model.rank = rank + 1 #setting post hoc
                model.best_model = best_model_type
                if rank != 0: #change in AIC relative to best model
                    model.delta_aic = abs(min_aic - model.AIC)
                    model.delta_bic = abs(min_bic - model.BIC)
                else:
                    model.delta_aic = 'NA'
                    model.delta_bic = 'NA'
                all_models.append(model)
    else:
        all_models = model_list
    return all_models
    
def compile_models(model_list, aux_cols = None, whole_curves = False, show_estimates=False, sortby=['Species', 'Model_name']):

    main_cols = ["Species", "Model_name", "Trait", "B0", "E", "T_pk", "E_D", "E_D_L", "Est_Tpk",
                 "Max_response", "Slope", "Intercept", "R_Squared", "AIC", "BIC", "Delta_AIC", "Delta_BIC", "Best_Type", 
                 "Rank", "Corrected", "Number_Of_Data_Points", "Points_Before_Peak", "Points_After_Peak", "Number_Of_Variables"]

    bootstrap_cols = ["B0_max", "B0_min", "E_max", "E_min", "Tpk_max", "Tpk_min", "Response_max",
                      "Response_min", "ED_max", "ED_min", "EDL_max", "EDL_min", "TH_max", "TH_min", 
                      "THL_max", "THL_min", "Slope_max", "Slope_min", "Intercept_max", "Intercept_min"]

    uncertainty_cols = ['Tpk_std', 'Response_std', 'E_std']

    whole_curve_cols = ['Temperature', 'Response', 'Original_Data']

    estimate_cols = ["B0_est", "E_est", "E_D_est", "E_D_L_est", "T_pk_est"]
                      
    aux_cols = aux_cols or [] #Safe way to ensure this list is always empty when the function is initialised
    any_bootstrapped = any([i.bootstrapped for i in model_list]) #do any models contain bootstrap data
    any_uncertainty = any([i.uncertainty_estimated for i in model_list]) #do any models contain uncertainty data

    #False * list = empty list so we can construct the header easily
    col_names = main_cols + (any_uncertainty * uncertainty_cols) + (show_estimates * estimate_cols) + \
                (any_bootstrapped * bootstrap_cols) + aux_cols + (whole_curves * whole_curve_cols) 

    rows = [] #place to store each row in the form (index, data)
    
    for model in model_list:
        model_details = [model.species_name, model.model_name, model.trait]
        model_est_params = get_estimated_parameters(model)
        
        fit_statistics =  [model.R2, model.AIC, model.BIC, model.delta_aic, model.delta_bic, 
                            model.best_model, model.rank, model.response_corrected, model.ndata,
                            getattr(model, 'points_before_peak', 'NA'), getattr(model, 'points_after_peak', 'NA'), model.nvarys]  
             
        model_output = model_details + model_est_params + fit_statistics

        if any_uncertainty:
            uncertainty_data = get_uncertainty_data(model)  
        else: 
            uncertainty_data = []

        if show_estimates:
            original_estimates = get_estimates(model) 
        else:
            original_estimates = []

        if any_bootstrapped:
            bootstrap_data = get_bootstrapped_parameters(model)  
        else: 
            bootstrap_data = []
        
        if whole_curves: #output the entire smooth growth curve
            #Output model curve
            for i, temp in np.ndenumerate(model.smooth_x):
                resp = model.smooth_y[i]
                key = model.index + '_model_point_' + str(i[0]) #needs to be unique
                
                entries = model_output + uncertainty_data + original_estimates + bootstrap_data + model.aux_parameters_values + [temp, resp, 'False']
                rows.append((key, entries))

            #Output original data points as well
            for i, temp in np.ndenumerate(model.temps):
                resp = model.responses[i]
                key = model.index + '_orig_point_' + str(i[0]) #needs to be unique
                
                entries = model_output + uncertainty_data + original_estimates + bootstrap_data + model.aux_parameters_values + [temp, resp, 'True']
                rows.append((key, entries))   
                
        else:
            entries = model_output + uncertainty_data + original_estimates + bootstrap_data + model.aux_parameters_values
                   
            row = (model.index, entries)
            rows.append(row)
    
    if len(rows[0]) != len(col_names):
        print('Warning: length of col names and rows is different, did you remember to include the additional columns when calling this function?')

    df = pd.DataFrame.from_items(rows,  orient='index', columns=col_names) #create the dataframe from a list of rows
    df = df.sort_values(sortby).fillna('NA') #replace any blanks with 'na'
    df = df.replace('#NAME?', 'NA')
 
    return df

def get_estimated_parameters(model):
    "Safely get model parameters, getattr will return 'NA' if parameters are not found"
    final_B0 = getattr(model, 'final_B0', "NA") #If attribute doesn't exist returns NA
    final_E = getattr(model, 'final_E', "NA")
    final_T_pk = getattr(model, 'final_T_pk', "NA")
    final_estimated_T_pk = getattr(model, 'tpk_est', "NA")
    final_max_response = getattr(model, 'max_response_est', "NA")
    final_E_D = getattr(model, 'final_E_D', "NA")
    final_E_D_L = getattr(model, 'final_E_D_L', "NA")
    final_T_H = getattr(model, 'final_T_H', "NA")
    final_T_H_L = getattr(model, 'final_T_H_L', "NA")
    final_slope = getattr(model, 'slope', "NA")
    final_intercept = getattr(model, 'intercept', "NA")

    model_parameters = [final_B0, final_E, final_T_pk, final_E_D, final_E_D_L, 
                        final_estimated_T_pk, final_max_response, final_slope, 
                        final_intercept]

    return model_parameters

def get_estimates(model):
    "Safely get model parameters, getattr will return 'NA' if parameters are not found"
    est_B0 = getattr(model, 'B0', "NA") #If attribute doesn't exist returns NA
    est_E = getattr(model, 'E_init', "NA")
    est_E_D = getattr(model, 'E_D_init', "NA")
    est_E_D_L = getattr(model, 'E_D_L_init', "NA")
    est_T_pk = getattr(model, 'T_pk', "NA")

    parameter_estimates = [est_B0, est_E, est_E_D, est_E_D_L, est_T_pk]

    return parameter_estimates

def get_bootstrapped_parameters(model):
    "Get the estimated parameters of a bootstrapped model safely"
    final_B0_max = getattr(model, 'final_B0_max', "NA")
    final_E_max = getattr(model, 'final_E_max', "NA")
    final_estimated_T_pk_max = getattr(model, 'tpk_est_max', "NA")
    final_max_response_max = getattr(model, 'max_response_est_max', "NA")
    final_E_D_max = getattr(model, 'final_E_D_max', "NA")
    final_E_D_L_max = getattr(model, 'final_E_D_L_max', "NA")
    final_T_H_max = getattr(model, 'final_T_H_max', "NA")
    final_T_H_L_max = getattr(model, 'final_T_H_L_max', "NA")
    final_slope_max = getattr(model, 'final_slope_max', "NA")
    final_intercept_max = getattr(model, 'final_intercept_max', "NA")

    final_B0_min = getattr(model, 'final_B0_min', "NA")
    final_E_min = getattr(model, 'final_E_min', "NA")
    final_estimated_T_pk_min = getattr(model, 'tpk_est_min', "NA")
    final_max_response_min = getattr(model, 'max_response_est_min', "NA")
    final_E_D_min = getattr(model, 'final_E_D_min', "NA")
    final_E_D_L_min = getattr(model, 'final_E_D_L_min', "NA")
    final_T_H_min = getattr(model, 'final_T_H_min', "NA")
    final_T_H_L_min = getattr(model, 'final_T_H_L_min', "NA")        
    final_slope_min = getattr(model, 'final_slope_min', "NA")      
    final_intercept_min = getattr(model, 'final_intercept_min', "NA")

    bootstrap_parameters = [final_B0_max, final_B0_min, final_E_max, final_E_min, final_estimated_T_pk_max,
                      final_estimated_T_pk_min, final_max_response_max, final_max_response_min, 
                      final_E_D_max, final_E_D_min, final_E_D_L_max, final_E_D_L_min, final_T_H_max,
                      final_T_H_min, final_T_H_L_max, final_T_H_L_min, final_slope_max, final_slope_min,
                      final_intercept_max, final_intercept_min]

    return bootstrap_parameters

def get_uncertainty_data(model):
    "Get the estimated uncertainty of a model safely"
    final_tpk_uncertainty = getattr(model, 'tpk_uncertainty', "NA")
    final_response_uncertainty = getattr(model, 'response_uncertainty', "NA")
    final_E_uncertainty = getattr(model, 'E_uncertainty', "NA")

    uncertainty_parameters = [final_tpk_uncertainty, final_response_uncertainty, final_E_uncertainty]

    return uncertainty_parameters

#Weighted arithmetic mean
def weighted_amean(data, weights, normalise_weights = False, geometric = False):
    data, weights =  clean_paired_data(data, weights)
    
    if np.all(np.isnan(data)) or np.all(np.isnan(weights)):
        return np.nan
    
    if normalise:
        weights = normalise(weights)
    
    a_mean = np.average(data, weights=weights)

    if geometric:
        return np.exp(a_mean)
    return a_mean    
    
#Weighted Std Dev
def weighted_std(data, weights, normalise_weights = False, geometric = False):
    data, weights =  clean_paired_data(data, weights)
    
    if len(data) == 1:
        return weights

    if np.all(np.isnan(data)) or  np.all(np.isnan(weights)):
        return np.nan    
    
    if normalise_weights:
        weights = normalise(weights)

    weights_sum = weights.sum()

    if geometric: 
        weighted_mean = np.exp(np.average(np.log(data), weights=weights))
    else:
        weighted_mean = np.average(data, weights=weights)

    wvar = (weights*(data-weighted_mean)**2).sum() / weights_sum 
    wsdev = np.sqrt(wvar) 
   
    if geometric:
        return np.exp(wsdev)
    return wsdev

#Normalise and array
def normalise(arr):
    #Convert uncertainties to weights
    arr = np.array(arr, dtype=np.float64)
    
    inverse = 1 / (arr + 1) #Add a small value so values with no uncertainty have highest weight
    finite_vals = inverse[np.isfinite(inverse)] #infinite values are largest but make poor divisorss
    
    if finite_vals.any():
        maximum_val = finite_vals.max()
        return inverse / maximum_val
    else:
        return 0

#Geometric mean robust to 0s
def gmean(data):
    data = np.array(data, dtype=np.float64)
    
    if np.all(np.isnan(data)):
        return np.nan
        
    orig_length = len(data)
    data = data[np.isfinite(data)]

    log_data = np.log(data)
    result = np.exp(sum(log_data / orig_length))
    
    return result
    
def clean_paired_data(arr1, arr2):
    "Remove nans from two arrays in the same positions"
    arr1, arr2 = np.array(arr1), np.array(arr2)

    #Convert the arrays into a singe x by 2 matrix and transpose into a 2 by x matrix
    paired = np.transpose([arr1, arr2])
    
    #Remove any rows containing an invalid entry
    cleaned = np.ma.compress_rows(np.ma.fix_invalid(paired))
    
    if np.any(cleaned):
        #Transpose back into x by 2 and index for the original arrays
        cleaned = np.transpose(cleaned)
        arr1_clean, arr2_clean = cleaned[0], cleaned[1]
    else:
        #If the array turns out to be totally invalid return nan
        return np.array([np.nan]), np.array([np.nan])
    
    return arr1_clean, arr2_clean

def randomise_estimates(parameters, N=100, scale_levels=1):
    "Randomise the starting parameters of a model fit to try and improve NLS fitting"
    parameters_list_randomised = []
    scale_length = int(np.ceil(N / scale_levels))

    #extract original estimates
    original_B0 = parameters.B0
    original_E_init = parameters.E_init
    original_E_D_init = parameters.E_D_init
    original_E_D_L_init = parameters.E_D_L_init
    original_T_H_L = parameters.T_H_L
    original_T_H = parameters.T_H
    original_T_pk = parameters.T_pk

    #rachet up the deviation as we go through to try and get a good fit
    for i in range(scale_levels):
        scale = i + 1
        #Sample from the normal distribution
        B0_vals = np.random.normal(original_B0,                 scale, size=scale_length)
        E_init_vals = np.random.normal(original_E_init,         scale, size=scale_length)
        E_D_init_vals = np.random.normal(original_E_D_init,     scale, size=scale_length)
        E_D_L_init_vals = np.random.normal(original_E_D_L_init, scale, size=scale_length)
        T_H_L_vals = np.random.normal(original_T_H_L,           scale, size=scale_length)
        T_H_vals = np.random.normal(original_T_H,               scale, size=scale_length)
        T_pk_vals = np.random.normal(original_T_pk,             scale, size=scale_length)

        parameters_list = [deepcopy(parameters) for i in range(scale_length)]

        #would do this as a list comprehension but its being weird and I'd rather watch the olympic tennis final right now
        for i, params in enumerate(parameters_list):
            params.set_estimates(B0=B0_vals[i], 
                                 E_init=E_init_vals[i],
                                 E_D_init = E_D_init_vals[i],
                                 E_D_L_init = E_D_L_init_vals[i],
                                 T_H_L = T_H_L_vals[i],
                                 T_H = T_H_vals[i],
                                 T_pk = T_pk_vals[i])

            parameters_list_randomised.append(params)

    return parameters_list_randomised

def find_linear_arrhenius(parameters, show = False):
    """find the most arrhenius like portion of a data distribution

    See Pawar et al 2016 for details
    """

    print('\t Attempting to linearise BA residuals')

    temps, responses = parameters.temps, parameters.responses
    temps, responses = clean_paired_data(temps, responses)

    #Create a temperature range to check over. Don't want to remove too much data!
    max_temp = max(temps)
    min_temp = max_temp - (max_temp - min(temps))
    peaktemp = temps[np.argmax(responses)]
    min_coefficient = None

    #Store the best parameters.
    best_parameters = deepcopy(parameters)

    while (max_temp - min_temp) > 10 and max_temp > peaktemp and len(temps) >= 8:
        #calc  arrhenius plot
        temps, responses = parameters.temps, parameters.responses

        #remove nans
        temps, responses = clean_paired_data(temps, responses)

        temps = 1 / temps
        responses = np.log(responses)

        #BA model again
        model = np.polyfit(temps, responses, 1)
        model_fn = np.poly1d(model)

        residuals = responses - (model_fn(temps))

        #Fit a 2nd order polynomial to residuals
        coefficient = np.polyfit(temps, residuals, 2)
        fitted = np.poly1d(coefficient)

        #Find the minimum possible X^2 coefficient
        if not min_coefficient or (abs(min_coefficient) > abs(coefficient[0]) and coefficient[1] < 0):

            if show:
                plt.plot(temps, responses, '.')
                plt.plot(temps, model_fn(temps))
                plt.plot(temps, fitted(temps), '-')
                plt.show()
                plt.close()

            print('\n\tFound new optima at {0:.2f}C'.format(max_temp))
            best_parameters = deepcopy(parameters)
            min_coefficient = abs(coefficient[0])
        
        #Find the temperature difference between current max to avoid pointless looping
        difference = max_temp - max(parameters.temps)

        #Reduce the temperature cutoff and try again
        max_temp -= difference

        #Remove weight, temperature and response values above the max temp.
        parameters.keep_between_temps(min_temp - 1, max_temp)

    return best_parameters

def plot_many_thermal(arr, path, legend = None, title_extra = ''):
    """semi intelligently plot several thermal model objects"""
    lines = itertools.cycle(('-', '--', '-.', ':'))

    f = plt.figure()
    sns.set_style("ticks", {'axes.grid': True})
    ax = f.add_subplot(111)
    
    species_names = list(set([i.species_name for i in arr]))
    title = ' '.join(species_names) + title_extra
    
    #Set a name for each line
    if not legend:
        legend = [i.model_name_short for i in arr]

    #Get raw data which models were fitted to
    data_x = [i.temps for i in arr]
    data_y = [i.responses for i in arr]

    #Scale the graph manually
    data_y_max = [i.max() for i in data_y]
    
    #Get smooth line data
    curve_x  = [i.smooth_x for i in arr]
    curve_y = [i.smooth_y for i in arr]

    #Get an e val for each line
    e_vals = [i.final_E for i in arr]
    
    weights_x, weights_y = [], []
    
    for i in arr:
        try: 
            weights_x.append(i.x_weights)
            weights_y.append(i.y_weights)
        except:
            pass
           
    ax.set_ylim((0, max(data_y_max) * 1.1))

    #Plot actual observed data as a scatterplot
    for x, y in zip(data_x, data_y):
        plt.plot(x, y, marker='o', ms=3, linestyle='None', color='green', alpha=0.7, zorder=1)
        
    #Plot fitted curve
    for x, y, name, e in zip(curve_x, curve_y, legend, e_vals):
        plt.plot(x, y, marker='None', linewidth=3, linestyle=next(lines), zorder=3, label = '{}, e = {:.2f}'.format(name, e))
    
    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    plt.title(title)

    plt.savefig(path + title_extra, bbox_inches='tight')
    plt.close()

def plot_weights_histograms(model, path, legend = None):

    title = model.species_name + '_' + model.model_name_short + '_weight_distribution'

    bins = np.linspace(0, 1, 20)
    weights_x, weights_y = model.x_weights, model.y_weights

    f = plt.figure()
    sns.set_style("ticks", {'axes.grid': True})
    ax = f.add_subplot(111)

    plt.hist(weights_x, bins, alpha=0.5, label='x')
    plt.hist(weights_y, bins, alpha=0.5, label='y')

    plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    plt.title(title)

    path = path + '_' + model.model_name_short + '_weight_distribution'

    plt.savefig(path, bbox_inches='tight')
    plt.close()
    
