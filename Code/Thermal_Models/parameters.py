import numpy as np
import pandas as pd
import os
import re

from scipy import stats, integrate

class estimate_parameters:
    """
    This class estimates all of metabolic parameters which are required as starting points for the least squared fitting of the models themselves.
    It also extracts useful data from the database which is passed to the models.
    
    I'd ideally like to split this into two parts in the future. One handling estimation, the other handling auxillary information.
    Really some of this could just be a named tuple.  
    """
    
    #Going to define a bunch of parameters here which can then be overwritten by passing flags to __init__
    
    k = 8.62e-5 #Boltzmann constant
    Tref = 273.15 #Reference temperature - 0C
    
    Trait_Col_Name = 'StandardisedTraitName' #trait name
    X_vals_col_name = 'ConTemp' #column name to pull x values from
    Y_vals_col_name = 'StandardisedTraitValue' #column name to pull y vals from
    uncertainty_x_col_name = None
    uncertainty_y_col_name = None

    x_val_conversion = 60 * 60 * 24
    
    species_name = ''
    
    full_strain_col_name = 'Consumer'
    genus_level_col_name = 'ConGenus'
    species_level_col_name = 'ConSpecies'
    
    species_data = True
    is_celcius = True #Is the input temps in celcius
    
    def __init__(self, data=None, aux_parameters_names = [] , flags = {}):
       
        for k, v in flags.items(): #flags will overwrite the values above allowing for flexible databases
            setattr(self, k, v)
       
        self.aux_parameters_names = aux_parameters_names
        
        if data is not None:   
            self.data = self.clean_dataset(data)
            self.get_ancillary_info()
            self.trait = self.get_single_val(self.Trait_Col_Name)
        
            self.temps = self.get_column('K') #Temperatures
            self.responses = self.get_column('Cor_Trait_Value') #Growth rates

            self.set_name()
            self.estimate_all()
    
    def estimate_all(self, boot=False):
        "Estimate all data points, this is outside __init__ so we can call it again when we bootstrap"
        self.get_T_pk() #Assign value for T - peak
        self.calc_slopes()
        self.estimate_E_init()
        self.estimate_low_temp_correction()
        self.estimate_high_temp_correction()
        self.estimate_B0()

        if not boot: #if we're bootstrapping don't reset the weights as this is iterative
            self.set_residual_weights()
        
    def get_single_val(self, item):
        "Return the top value of a column in a pandas dataframe"
        return self.data[item][0]
        
    def get_column(self, item, remove_nan = True):
        "Safe way to get a numerical column from a pandas dataframe"
        col = self.data[item]
        col = col.replace('NA', np.nan) #replace with nans so we can do maths without worring about strings
        col = col.as_matrix() #convert to numpy array
        if remove_nan:
            return col[~np.isnan(col)] #Remove Nan values
        return col
        
    def resample_data(self):
        "resample so we can bootstrap"
        bootstrap_N = len(self.temps)
        
        indices = np.random.choice(bootstrap_N, bootstrap_N) #create a vector of list indexes the same length as the original list
        
        self.temps = self.temps[indices] #Resample keeping x and y coupled
        self.responses = self.responses[indices]

        if isinstance(self.uncertainty_x, (np.ndarray, np.generic)):
            self.uncertainty_x = self.uncertainty_x[indices]

        if isinstance(self.uncertainty_y, (np.ndarray, np.generic)):
            self.uncertainty_y = self.uncertainty_y[indices]

    def set_estimates(self, **kwargs):
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

    def keep_between_temps(self, low_temp, high_temp):
        "Remove all TPK and UPK pairs not within a certain temperature"
        indices = np.where(np.logical_and(self.temps > low_temp, self.temps < high_temp)) #create a truth array of values within the desired range

        self.temps = self.temps[indices] #Resample keeping x and y coupled
        self.responses = self.responses[indices]

        if isinstance(self.uncertainty_x, (np.ndarray, np.generic)):
            self.uncertainty_x = self.uncertainty_x[indices]

        if isinstance(self.uncertainty_y, (np.ndarray, np.generic)):
            self.uncertainty_y = self.uncertainty_y[indices]

    def clean_dataset(self, data):
        "Normalise each dataset"
        #Transform temps to kelvin
        if self.is_celcius:
            data['K'] = data[self.X_vals_col_name] + 273.15
        else:
            data['K'] = data[self.X_vals_col_name]
        
        # Convert corrected value from s^-1 to d^-1
        data['Cor_Trait_Value'] = data[self.Y_vals_col_name] * self.x_val_conversion # Convert corrected value from s^-1 to d^-1
        
        #If any trait values are negative then subtract the smallest value to normalise
        minimum_temp_value  = data['K'].min()
        minimum_trait_value  = data['Cor_Trait_Value'].min()
        
        if minimum_trait_value <= 0:
            data['Cor_Trait_Value'] -= minimum_trait_value - 10E-10 #Get rid of any 0s
        
        if minimum_temp_value <= 0:
            data['K'] -= minimum_temp_value - 10E-10 #Get rid of any 0s
            
        return data
        
    def get_ancillary_info(self):
        "Get information on each curve to include in the summary"
        self.aux_parameters_values = [self.data[aux_parameter][0] for aux_parameter in self.aux_parameters_names]  
        
    def get_T_pk(self):
        "Find the temperature at which maximum response is observed"
        self.Tpk_row = self.responses.argmax() #Index of max response
        self.T_pk = self.temps[self.Tpk_row] #Temperature at which rate is maximum
        
    def calc_slopes(self):
        "Slice the data to find the upwards and downwards slopes in the dataset"
        self.upslope_x, self.downslope_x = self.temps[:self.Tpk_row + 1], self.temps[self.Tpk_row:]
        self.upslope_y, self.downslope_y = self.responses[:self.Tpk_row + 1], self.responses[self.Tpk_row:]  

    def estimate_E_init(self):
        "Estimate energy value using the slope of the values to the peak of an arrhenius plot"
        if len(self.upslope_x) > 1:
            x = 1 / (self.k * self.upslope_x)
            y = np.log(self.upslope_y) 

            try:
                slope,*vals = stats.theilslopes(x,y, 0.9) #maybe more robust given noisy data? 
            except:
                slope,*vals = stats.linregress(x,y) 

            self.E_init = abs(slope)
        else:
            self.E_init = 0.6 #Default value

    def estimate_high_temp_correction(self):
        "Estimate energy value using the slope of the values to the peak of an arrhenius plot"
        
        if len(self.downslope_x) > 1:
            #Estimate ED
            x = 1 / (self.k * self.downslope_x)
            y = np.log(self.downslope_y)
            
            try:
                slope, *vals = stats.theilslopes(x,y, 0.9) #maybe more robust given noisy data? 
            except:
                slope, *vals = stats.linregress(x,y)

            self.E_D_init = slope + self.E_init

            #Estimate TH
            downslope_diff_x = self.downslope_x[1:]
            downslope_diff_y = np.diff(self.downslope_y)
            
            max_change_index = np.argmin(downslope_diff_y)
            self.T_H = self.T_pk + ((downslope_diff_x[max_change_index] - self.T_pk) / 2)
        else: 
            self.E_D_init = self.E_init * (4) #Default value
            self.T_H = self.T_pk + 3

    def estimate_low_temp_correction(self):
        "Estimate energy value using the slope of the values to the peak of an arrhenius plot"
        if len(self.upslope_x) > 3:
            #Estimste THL
            upslope_diff_x = self.upslope_x[:-2]
            upslope_diff_y = np.diff(self.upslope_y[:-1])

            max_diff_index = np.argmax(upslope_diff_y)
            min_diff_index = np.argmin(upslope_diff_y)
            self.T_H_L = (upslope_diff_x[max_diff_index] + upslope_diff_x[min_diff_index]) / 2
        else: 
            self.T_H_L = self.T_pk - 10

        if len(self.upslope_x) > 5:
            #estimate EDL
            x = np.array_split(self.upslope_x, 3)[0]
            y = np.array_split(self.upslope_y, 3)[0]

            ahrr_x = 1 / (self.k * x)
            ahrr_y = np.log(y)

            
            try: 
                slope, *vals = stats.theilslopes(ahrr_x, ahrr_y, 0.9) #maybe more robust given noisy data? 
            except:
                slope, *vals = stats.linregress(ahrr_x, ahrr_y)

            self.E_D_L_init = slope + self.E_init
        else: 
            self.E_D_L_init = self.E_init * (-2) #Default value

    def estimate_B0(self):
        "Returns the response at the tempetature closest to Tref"
        closest_T_index = abs(self.temps - (self.Tref + 25)).argmin() #index of data point closest to 25C
        closest_T_response = self.responses[closest_T_index]
        if closest_T_response != 0: 
            self.B0 = np.log(closest_T_response)
        else:
            self.B0 = np.log(self.responses.max()) 
    
    def set_residual_weights(self):
        if self.uncertainty_x_col_name: #set residual weights, other 1 = no weighting
            self.uncertainty_x = self.get_column(self.uncertainty_x_col_name, remove_nan = False)
        else:
            self.uncertainty_x = 1

        if self.uncertainty_y_col_name: #set residual weights, other 1 = no weighting
            self.uncertainty_y = self.get_column(self.uncertainty_y_col_name, remove_nan = False)
        else:
            self.uncertainty_y = 1
                            
    def set_name(self):
        "Set species name to be applied to plot title"
        if self.species_name == '' and isinstance(self.species_name, str):
            genus = self.get_single_val(self.genus_level_col_name)
            species = self.get_single_val(self.species_level_col_name)
            consumer = self.get_single_val(self.full_strain_col_name)
        
            #Use this to remove pseudoreplicates
            if pd.isnull(species) or not self.species_data:
                self.species_name = consumer #if no species is available we have to use consumer
            else:
                self.species_name = ' '.join([genus, species])
            try:    
                self.species_name = self.species_name[0].upper() + self.species_name[1:].lower() #Usually a genus and species name so this should be correct in most cases
            except TypeError:
                print('Warning, no name found at this level for group')

    def __deepcopy__(self, memodict={}):
        """Deep copy is very slow so I'm defining a custom method to speed it up a bit!"""
        copy_object = estimate_parameters()

        copy_object.temps = self.temps
        copy_object.responses = self.responses
        copy_object.species_name = self.species_name
        copy_object.trait = self.trait
        copy_object.B0 = self.B0
        copy_object.E_init = self.E_init
        copy_object.E_D_init = self.E_D_init
        copy_object.E_D_L_init = self.E_D_L_init
        copy_object.T_H_L = self.T_H_L
        copy_object.T_H = self.T_H
        copy_object.T_pk = self.T_pk
        copy_object.uncertainty_x = self.uncertainty_x
        copy_object.uncertainty_y = self.uncertainty_y
        copy_object.aux_parameters_names = self.aux_parameters_names
        copy_object.aux_parameters_values = self.aux_parameters_values

        copy_object.Trait_Col_Name = self.Trait_Col_Name
        copy_object.X_vals_col_name = self.X_vals_col_name
        copy_object.Y_vals_col_name = self.Y_vals_col_name
        copy_object.uncertainty_x_col_name = self.uncertainty_x_col_name
        copy_object.uncertainty_y_col_name = self.uncertainty_y_col_name
        copy_object.x_val_conversion = self.x_val_conversion
        copy_object.species_name = self.species_name
        copy_object.full_strain_col_name = self.full_strain_col_name
        copy_object.genus_level_col_name = self.genus_level_col_name
        copy_object.species_level_col_name = self.species_level_col_name
        copy_object.species_data = self.species_data
        copy_object.is_celcius = self.is_celcius

        return copy_object
    
    def __str__(self):
        vars = [self.species_name, self.temps, self.responses, self.trait, self.B0,
                self.E_init, self.E_D_init, self.E_D_L_init, self.T_H_L, self.T_H, self.T_pk, self.Tpk_row]
        
        text = """
        ----------------------
        {0[0]}
        Trait: {0[3]}
        
        Estimates: 
        B0: {0[4]:.2f}
        E: {0[5]:.2f}
        ED: {0[6]:.2f}
        EDL: {0[7]:.2f}
        THL: {0[8]:.2f}
        TH: {0[9]:.2f}
        TPK: {0[10]:.2f}
        TPK Row: {0[11]:.2f}
        """.format(vars)
        
        return text     