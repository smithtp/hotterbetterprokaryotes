import numpy as np

from lmfit import Parameters
from .base_model import physiological_growth_model #All classes in this file inherit methods from physiological_growth_model

class schoolfield_two_factor_tpk(physiological_growth_model):
    "Schoolfield model using T_pk as a substitute for T_H"
    model_name = "schoolfield_two_factor_tpk"
    model_name_short = "SCH_TF_TPK"
        
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_TF_TKP" #ID for the model
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.get_final_values()
        self.assess_model()
        self.get_stderrs()
        self.find_peak()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                         Name,       Start,           Can_Vary,  Lower,           Upper
        self.parameters.add_many(('B0_start', self.B0,         True,      -np.inf,         np.inf,       None),
                                 ('E',        self.E_init,     True,       10E-3,          np.inf,       None),
                                 ('E_D',      self.E_D_init,   True,       10E-10,         np.inf,       None),
                                 ('T_pk',     self.T_pk,       True,       273.15,         np.inf,       None))


    def fit(self, parameter_vals, temps):
        "Fit a schoolfield curve to a list of temperature values"
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        T_pk = parameter_vals['T_pk'] #Temperature at which peak response is observed  
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref))) /\
                         (1 + (E/(E_D - E)) * np.exp(E_D / self.k * (1 / T_pk - (1 / temps)))))

        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_pk = values['T_pk']   

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_T_pk_stderr = self.model.params['T_pk'].stderr
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init,
                self.final_E, self.T_pk, self.final_T_pk, self.E_D_init, 
                self.final_E_D, self.R2, self.AIC, self.BIC]
                
        text = """
        ---Schoolfield Two Factor Model With Explicit TPK---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}
        
        T Peak est = {0[5]:.2f}
        T Peak final =  {0[6]:.2f}
        
        E_D est = {0[7]:.2f}
        E_D final = {0[8]:.2f}
        
        R2: = {0[9]:.2f}
        AIC = {0[10]:.2f}
        BIC = {0[11]:.2f}
        
        """.format(vars)
        return text

class schoolfield_two_factor(physiological_growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # This seems to generate pretty much identical results to the two factor model, but throws its toys out the pram sometimes.
     # In this version T_H is the temperature at which half the enzyme is denatured by heat stress
    model_name = "schoolfield_two_factor"
    model_name_short = "SCH_TF"
    
    def fit_from_parameters(self, est_parameters, index, fast=False):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_TF" #Used to name plot graphics file
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.get_final_values()
        self.assess_model()
        self.get_stderrs()
        self.find_peak()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                         Name,       Start,           Can_Vary,  Lower,           Upper
        self.parameters.add_many(('B0_start', self.B0,         True,     -np.inf,          np.inf,      None),
                                 ('E',        self.E_init,     True,      0.05,            np.inf,      None),
                                 ('E_D',      self.E_D_init,   True,      0,               np.inf,      None),
                                 ('T_H',      self.T_H,        True,      self.T_pk,       273.15+170,  None))

    def fit(self, parameter_vals, temps):
        "Fit a schoolfield curve to a list of temperature values"
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        T_H = parameter_vals['T_H'] #Temperature at which half od enzzymes are denatured
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))\
                        /(1 + np.exp((E_D / self.k) * (1 / T_H - (1 / temps)))))
        
        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_H = values['T_H']   
        
    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_T_stderr = self.model.params['T_H'].stderr
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, self.T_pk, self.T_H, 
                self.final_T_H, self.E_D_init, self.final_E_D, self.R2, self.AIC, self.BIC]
        text = """
        ---Schoolfield Two Factor Model---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}
        
        TPK = {0[5]:.2f}
        T H est = {0[6]:.2f}
        T H final =  {0[7]:.2f}
        
        E_D est = {0[8]:.2f}
        E_D final = {0[9]:.2f}
        
        R2: = {0[10]:.2f}
        AIC = {0[11]:.2f}
        BIC = {0[12]:.2f}
     
        """.format(vars)
        return text        

class schoolfield_full(physiological_growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # In this version T_H_L is a low temperature enzyme inactivation constant (as if having a high temp one wasn't fun enough already)
    model_name = "schoolfield_full"
    model_name_short = "SCH_F"
    
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_F" #Used to name plot graphics file
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.get_final_values()
        self.assess_model()
        self.get_stderrs()
        self.find_peak()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression, note additional T_H_L parameter"
        self.parameters = Parameters()
        
        #                         Name,       Start,                Can_Vary,  Lower,           Upper
        self.parameters.add_many(('B0_start', self.B0,               True,   -np.inf,         np.inf,           None),
                                 ('E',        self.E_init,           True,   0.05,            np.inf,           None), 
                                 ('E_D',      self.E_D_init,         True,   10E-10,          np.inf,           None),
                                 ('E_D_L',    self.E_D_L_init,       True,   -np.inf,         -10E-10,          None),
                                 ('T_H',      self.T_H,              True,   273.15-70,       273.15+170,       None),
                                 ('T_H_L',    self.T_H_L,            True,   273.15-70,       273.15+170,       None))

    def fit(self, parameter_vals, temps):
        "Fit a schoolfield curve to a list of temperature values"
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        E_D_L = parameter_vals['E_D_L'] #Energy of cold inactivation
        T_H = parameter_vals['T_H'] #Temperature at which half of enzymes are denatured
        T_H_L = parameter_vals['T_H_L'] #Temperature at which half of enzymes are cold inactivated
        
        #Put these here so we don't need to limit the parameters themselves - more flexible
        #TH must be greater than THL
        if T_H < (T_H_L + 1):
            T_H = T_H_L + 1

        #And conversely THL must be less than TH 
        if T_H_L > T_H - 1:
            T_H_L = T_H - 1

        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))    \
                        /(1 + np.exp((E_D_L / self.k) * (1 / T_H_L - (1 / temps))) +   \
                              np.exp((E_D / self.k) * (1 / T_H - (1 / temps))))
                          )
        
        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_E_D_L = values['E_D_L']
        self.final_T_H = values['T_H']
        self.final_T_H_L = values['T_H_L']

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_E_D_L_stderr = self.model.params['E_D_L'].stderr
        self.final_T_stderr = self.model.params['T_H'].stderr        
        self.final_T_H_L_stderr = self.model.params['T_H_L'].stderr    
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, 
                self.E_D_init , self.final_E_D, self.E_D_L_init, self.final_E_D_L,
                self.T_pk, self.T_H, self.final_T_H, self.T_H_L, self.final_T_H_L, 
                self.R2, self.AIC, self.BIC]
                
        text = """
        ---Schoolfield Three Factor Model---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}     
        
        E D est = {0[5]:.2f}
        E D final = {0[6]:.2f}
        
        E D L est = {0[7]:.2f}
        E D L final = {0[8]:.2f}
        
        TPK = {0[9]:.2f}
        T H est = {0[10]:.2f}
        T H final =  {0[11]:.2f}
        T H L est = {0[12]:.2f}
        T H L final =  {0[13]:.2f}  
        
        R2: = {0[14]:.2f}
        AIC = {0[15]:.2f}
        BIC = {0[16]:.2f}
     
        """.format(vars)
        return text                

class schoolfield_full_tpk(physiological_growth_model):
     # Schoolfield model using T_pk as a substitute for T_H, but including the full low temp inactivation terms
     # In this version T_H_L is a low temperature enzyme inactivation constant (as if having a high temp one wasn't fun enough already)
     # Behaves mostly the same as the original full schoolfield model, but there are instances where it behaves more like the simple formulation with explicit Tpk
    model_name = "schoolfield_full_tpk"
    model_name_short = "SCH_F_TPK"
    
    def fit_from_parameters(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_F_TPK" #Used to name plot graphics file
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.get_final_values()
        self.assess_model()
        self.get_stderrs()
        self.find_peak()
  
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression, note additional T_H_L parameter"
        self.parameters = Parameters()
        
        #                         Name,       Start,                Can_Vary,  Lower,           Upper
        self.parameters.add_many(('B0_start', self.B0,               True,     -np.inf,         np.inf,            None),
                                 ('E',        self.E_init,           True,     10E-3,           np.inf,            None), #E lower to 0.05 to stop it tending to 0
                                 ('E_D',      self.E_D_init,         True,     10E-10,          np.inf,            None),
                                 ('E_D_L',    self.E_D_L_init,       True,     -np.inf,         -10E-10,           None),
                                 ('T_H_L',    self.T_H_L,            True,     273.15-70,       np.inf,            None),
                                 ('T_pk',     self.T_pk,             True,     273.15,          np.inf,            None))


    def fit(self, parameter_vals, temps):
        "Fit a schoolfield curve to a list of temperature values"
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        E_D_L = parameter_vals['E_D_L'] #Energy of cold inactivation
        T_H_L = parameter_vals['T_H_L'] #Temperature at which half of enzymes are cold inactivated
        T_pk = parameter_vals['T_pk'] #Temperature at which peak response is observed 
        
        #Put these here so we don't need to limit the parameters themselves - more flexible
        #THL must be less than TPK
        if T_H_L >  T_pk - 1:
            T_H_L =  T_pk - 1
        
        #TPK must be more than THL
        if T_pk < T_H_L + 1:
            T_pk = T_H_L + 1

        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref))) / \
                         (1 + np.exp((E_D_L / self.k) * (1 / T_H_L - (1 / temps))) +  \
                         (E/(E_D - E)) * np.exp(E_D / self.k * (1 / T_pk - (1 / temps)))))
        
        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_E_D_L = values['E_D_L']
        self.final_T_H_L = values['T_H_L']
        self.final_T_pk = values['T_pk']   

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_E_D_L_stderr = self.model.params['E_D_L'].stderr        
        self.final_T_H_L_stderr = self.model.params['T_H_L'].stderr
        self.final_T_pk_stderr = self.model.params['T_pk'].stderr            
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, 
                self.E_D_init, self.final_E_D, self.E_D_L_init, self.final_E_D_L,
                self.T_pk, self.final_T_pk, self.T_H_L, self.final_T_H_L, 
                self.R2, self.AIC, self.BIC]
                
        text = """
        ---Schoolfield Three Factor Model with Tpk as an explicit parameter---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}     
        
        E D est = {0[5]:.2f}
        E D final = {0[6]:.2f}
        
        E D L est = {0[7]:.2f}
        E D L final = {0[8]:.2f}
        
        T Peak est = {0[9]:.2f}
        T Peak final = {0[10]:.2f}

        T H L est = {0[11]:.2f}
        T H L final =  {0[12]:.2f}  
        
        R2: = {0[13]:.2f}
        AIC = {0[14]:.2f}
        BIC = {0[15]:.2f}
     
        """.format(vars)
        return text                