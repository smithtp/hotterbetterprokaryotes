import numpy as np


from lmfit import Parameters, minimize
from scipy import stats
from .base_model import physiological_growth_model #All classes in this file inherit methods from physiological_growth_model

class boltzmann_arrhenius_two_weights(physiological_growth_model):
    """
    Boltzmann Arrhenius model with weights on both the X and Y residuals
    Quite a lot of rewriting here - not very standardised.
    """
    model_name = "boltzmann_arrhenius_two_weights"
    model_name_short = "BA_TW"
        
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.setup()
        self.re_estimate_e()
        self.index = str(index) + "_" + self.model_name_short #ID for the model
        self.set_parameters()
        self.fit_model()
        self.get_final_values()
        self.smooth()
        self.assess_model()
        self.get_stderrs()

    def setup(self):
        "Convert data to ahhrenius parameters, create a list of weights for each datapoint"
        self.inv_temps = 1 / (self.temps * self.k)
        self.log_responses = np.log(self.responses)

        self.x_weights = self.normalise_uncertainty(self.uncertainty_x) #set the nans to 0 to avoid weirdness when these are summed
        self.y_weights = self.normalise_uncertainty(self.uncertainty_y) 
    
    def re_estimate_e(self):
        "As this is weird weighted regression weight using a normal regression" 
        E, B0, *args = stats.linregress(self.inv_temps , self.log_responses)
        
        self.E_init = -E * self.k
        self.B0 = B0

    def set_parameters(self):
        self.parameters = Parameters()
        #                         Name,       Start,           Can_Vary,  Lower,           Upper
        self.parameters.add_many(('B0_start', self.B0,         True,      -np.inf,         np.inf,       None),
                                 ('E',        -self.E_init,    True,      -np.inf,         np.inf,       None))

    def fit_model(self, method="leastsq"):
        "Least squares regression to minimise fit"
        self.model = minimize(self.get_residuals, 
                              self.parameters, 
                              args=(self.inv_temps, self.log_responses, self.x_weights, self.y_weights),
                              method=method)
                              
        self.R2 = 1 - np.var(self.model.residual) / np.var(self.responses)
        self.ndata = self.model.ndata
        self.nvarys = self.model.nvarys
   
    def get_residuals(self, params, temps, responses, weights_x, weights_y):
        "Called by fit model only, generates residuals using test values"
        parameter_values = params.valuesdict()
        x_residuals = self.fit_x(parameter_values, responses) - temps
        y_residuals = self.fit_y(parameter_values, temps) - responses

        x_residuals = weights_x * np.power(x_residuals, 2)
        y_residuals = weights_y * np.power(y_residuals, 2)

        weighted_residuals = np.sqrt(x_residuals + y_residuals)

        return weighted_residuals

    def fit_y(self, parameter_vals, temps):
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes

        fit = (E * temps) + B0

        return fit
        
    def fit_x(self, parameter_vals, responses):
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes

        fit = (responses - B0) / E
        
        return fit
        
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1) #Extrapolate a little 
        self.smooth_y = np.exp(self.final_B0) * np.exp(-self.final_E / (self.smooth_x * self.k)) 

    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = -values['E']

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init,
                self.final_E, self.R2, self.AIC, self.BIC]
                
        text = """
        ---Schoolfield Two Factor Model With Explicit TPK---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}

        R2: = {0[5]:.2f}
        AIC = {0[6]:.2f}
        BIC = {0[7]:.2f}
        
        """.format(vars)
        return text

#These are temporary until we work out which is better

class boltzmann_arrhenius_two_weights_non_linear(boltzmann_arrhenius_two_weights):
    model_name = "boltzmann_arrhenius_two_weights_non_linear"
    model_name_short = "BA_TWNL"

    def normalise_uncertainty(self, arr):
        """use a linear mapping of uncertainty to weights"""
        inverse = 1 / (arr + 1) #add a small number here so we don't divide by 0.
        finite_vals = inverse[np.isfinite(inverse)] #infinite values are largest but make poor divisors

        if finite_vals.any():
            maximum_val = finite_vals.max()
            weights = inverse / maximum_val
        else:
            weights = inverse

        weights[np.isnan(weights)] = 0
       
        return weights

class boltzmann_arrhenius_two_weights_linear(boltzmann_arrhenius_two_weights):
    model_name = "boltzmann_arrhenius_two_weights_linear"
    model_name_short = "BA_TWL"

    def normalise_uncertainty(self, arr):
        """use a non linear mapping of uncertainty to weights
        
        ISSUE: Weights lowest certainty point as 0
        """
        finite_vals = arr[np.isfinite(arr)] #infinite values are largest but make poor divisorss

        if finite_vals.any():
            maximum_val = finite_vals.max()

            normalised = arr / maximum_val

            weights = 1 - normalised
        else: 
            weights = arr

        weights[np.isnan(weights)] = 0

        return weights

class boltzmann_arrhenius_two_weights_normalise_only(boltzmann_arrhenius_two_weights):
    model_name = "boltzmann_arrhenius_two_weights_normalise_only"
    model_name_short = "BA_TWNO"

    def normalise_uncertainty(self, arr):
        """use a non linear mapping of uncertainty to weights
        
        ISSUE: Weights lowest certainty point as 0
        """
        finite_vals = arr[np.isfinite(arr)] #infinite values are largest but make poor divisorss

        if finite_vals.any():
            maximum_val = finite_vals.max()

            weights = arr / maximum_val

        else: 
            weights = arr

        weights[np.isnan(weights)] = 0

        return weights

class boltzmann_arrhenius_one_weight_linear(boltzmann_arrhenius_two_weights_normalise_only):
    model_name = "boltzmann_arrhenius_one_weight_linear"
    model_name_short = "BA_OWL"

    def get_residuals(self, params, temps, responses, weights_x, weights_y):
        "Called by fit model only, generates residuals using test values"
        parameter_values = params.valuesdict()

        y_residuals = self.fit_y(parameter_values, temps) - responses
        y_residuals = weights_y * np.power(y_residuals, 2)

        weighted_residuals = np.sqrt(y_residuals)

        return weighted_residuals

class boltzmann_arrhenius_one_weight_non_linear(boltzmann_arrhenius_two_weights_non_linear):
    model_name = "boltzmann_arrhenius_one_weight_non_linear"
    model_name_short = "BA_OWNL"

    def get_residuals(self, params, temps, responses, weights_x, weights_y):
        "Called by fit model only, generates residuals using test values"
        parameter_values = params.valuesdict()

        y_residuals = self.fit_y(parameter_values, temps) - responses
        y_residuals = weights_y * np.power(y_residuals, 2)

        weighted_residuals = np.sqrt(y_residuals)

        return weighted_residuals