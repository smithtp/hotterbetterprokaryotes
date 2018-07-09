import numpy as np

from lmfit import Parameters
from scipy import stats
from .base_model import physiological_growth_model #All classes in this file inherit methods from physiological_growth_model

class boltzmann_arrhenius(physiological_growth_model):
    model_name = "boltzmann_arrhenius"
    model_name_short = "BA"
    
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_BA" #ID for the model
        self.fit_model()
        self.smooth()
        self.assess_model()
        
    def fit_model(self): #Note this overwrites the method in the parent class so we totally cut the NLS component from this model
        x = 1 / (self.temps)
        y = np.log(self.responses)
        
        E, B0, r, p_value, std_err = stats.linregress(x, y)
        
        self.final_E = -E * self.k
        self.final_B0 = B0
        
        self.R2 = r * r
        
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1) #Extrapolate a little 
        self.smooth_y = np.exp(self.final_B0) * np.exp(-self.final_E / (self.smooth_x * self.k))        
        
    def assess_model(self):
        k = 2 #Number of variables
        n = len(self.temps) #Number of data points
        
        rss = self.R2 #Residual sum of squares
        
        self.ndata = n
        self.nvarys = 2
        self.AIC = 2 * k + n * np.log(rss / n)
        self.BIC = np.log(n)*(k + 1) + n * np.log(rss / n)
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.final_E, self.final_B0, self.R2, self.AIC, self.BIC]
                
        text = """
        --- Boltzmann Arrhenius ---
        {0[0]}
        
        E = {0[1]:.2f}
        B0 = {0[2]:.2f}
        
        R2: = {0[3]:.2f}
        AIC = {0[4]:.2f}
        BIC = {0[5]:.2f}
        
        """.format(vars)
        return text                
        
class lm(physiological_growth_model):
    "linear_model"
    model_name = "lm"
    model_name_short = "LM"
    
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_LM" #ID for the model
        self.fit_model()
        self.smooth()
        self.assess_model()
        
    def fit_model(self):
        self.slope, self.intercept, r, p_value, std_err = stats.linregress(self.temps, self.responses)
        
        self.R2 = r * r
        
        self.ndata = len(self.temps)
        self.nvarys = 2
        
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1) #Extrapolate a little 
        self.smooth_y = self.intercept + (self.smooth_x * self.slope)        
        
    def assess_model(self):
        k = 2 #Number of variables
        n = len(self.temps) #Number of data points
        
        rss = self.R2 #Residual sum of squares
        
        self.AIC = 2 * k + n * np.log(rss / n)
        self.BIC = np.log(n)*(k + 1) + n * np.log(rss / n)
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.slope, self.intercept, self.R2, self.AIC, self.BIC]
                
        text = """
        ---Linear Model (yawn)---
        {0[0]}
        
        Slope = {0[1]:.2f}
        Intercept = {0[2]:.2f}
        
        R2: = {0[3]:.2f}
        AIC = {0[4]:.2f}
        BIC = {0[5]:.2f}
        
        """.format(vars)
        return text 