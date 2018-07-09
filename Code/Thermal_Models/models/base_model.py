import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import re

from lmfit import Parameters, minimize, Minimizer, fit_report, conf_interval, printfuncs

class physiological_growth_model:
    k = 8.62e-5 #Boltzmann constant
    Tref = 273.15 #Reference temperature - 0C
    
    response_corrected = False #In cases where a peak in the curve tends to infinity this will be set to true. 
    bootstrapped = False
    uncertainty_estimated = False
    rank = 'NA' #When fitting multiple models to the same data we will use this to rank models
    
    def __init__(self):
        self.temps = np.array([])
        self.responses = np.array([])
        self.species_name = None
        self.trait = None
        self.B0 = None
        self.E_init = None
        self.E_D_init = None
        self.E_D_L_init = None
        self.T_H_L = None
        self.T_H = None
        self.T_pk = None
        self.aux_parameters_names = []
        self.aux_parameters_values = []
        
    def extract_parameters(self, est_parameters):
        self.temps = est_parameters.temps
        self.responses = est_parameters.responses
        self.species_name = est_parameters.species_name
        self.trait = est_parameters.trait
        self.B0 = est_parameters.B0
        self.E_init = est_parameters.E_init
        self.E_D_init = est_parameters.E_D_init
        self.E_D_L_init = est_parameters.E_D_L_init
        self.T_H_L = est_parameters.T_H_L
        self.T_H = est_parameters.T_H
        self.T_pk = est_parameters.T_pk
        self.uncertainty_x = est_parameters.uncertainty_x
        self.uncertainty_y = est_parameters.uncertainty_y
        self.aux_parameters_names = est_parameters.aux_parameters_names
        self.aux_parameters_values = est_parameters.aux_parameters_values
         
    def get_residuals(self, params, temps, responses):
        "Called by fit model only, generates residuals using test values"
        parameter_values = params.valuesdict()
        residuals = np.exp(self.fit(parameter_values, temps)) - responses

        return residuals
        
    def fit_model(self, method="leastsq"):
        "Least squares regression to minimise fit"
        self.mini = Minimizer(self.get_residuals, self.parameters, fcn_args=(self.temps, self.responses))

        self.model = self.mini.minimize(method='leastsq')

        """
        If all you care about is the fit use these lines in place of the last one, its slower though!
        start = self.mini.minimize(method='Nelder') #start with a slower but more accurate nelder method
        self.model = self.mini.minimize(method='leastsq', params=start.params) #use leastsq to generate standarderrors
        """

        self.R2 = 1 - np.var(self.model.residual) / np.var(self.responses)
        self.ndata = self.model.ndata
        self.nvarys = self.model.nvarys
        
    def assess_model(self):
        """Calculate the Akaike Information Criterion and Bayesian Information Criterion, using:

        - n:   number of observations
        - k:   number of parameters
        - rss: residual sum of squares
        """
        k = self.model.nvarys #Number of variables
        n = self.model.ndata #Number of data points
        rss = sum(np.power(self.model.residual, 2)) #Residual sum of squares

        #temp penalisation - would rather have models with covar
        self.AIC = 2 * k + n * np.log(rss / n)
        self.BIC = np.log(n)*(k + 1) + n * np.log(rss / n)

        if not (self.model.covar is not None and np.all(np.isfinite(self.model.covar))):
            self.AIC += abs(self.AIC)

        #temp penalisation - would rather have models without minimal E
        if self.final_E == 10E-3:
            self.AIC +=  0.2 * abs(self.AIC)

        #find data points before peak and after peak
        self.points_before_peak = len(self.temps[self.temps <= getattr(self, 'final_T_pk', 0)])
        self.points_after_peak = len(self.temps[self.temps > getattr(self, 'final_T_pk', 0)])
  
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1) #Extrapolate a little 
        self.smooth_y = np.exp(self.fit(self.model.params, self.smooth_x))
        
    def find_peak(self):
        """
        set the peak temperature and peak response of the main curve
        """ 
        temps = self.temps
        parameters = self.model.params

        self.tpk_est, self.max_response_est = self.differentiate_schoolfield(temps, parameters)
    
    def differentiate_schoolfield(self, temps, parameters):
        "Differentiate a schoolfield function numerically (actually interpolation for the moment)"
        tpk_est, max_response_est = 'NA', 'NA' #Set intitial values

        #if we estimated tpk directly use that
        if 'T_pk' in parameters and False: #temporraily disabled, fall back on interpolation
            if type(parameters) != dict:
                parameters = parameters.valuesdict()

            tpk_est = parameters['T_pk']
            x_vals = np.array([tpk_est])
            fit = np.exp(self.fit(parameters, x_vals))
            max_response_est = fit[0]

        #Otherwise use interpolation
        else:
            x_vals = np.arange(temps.min() - 5, temps.max() + 5, 0.1) #Extrapolate
            fit = self.fit(parameters, x_vals)   #Fit the model using the parameters

            peak_check_x = x_vals[np.isfinite(fit)] #remove the nans and infs so we can take the exponential
            peak_check_y = np.exp(fit[np.isfinite(fit)])

            #differences between values
            differential = np.diff(peak_check_y)

            turning_points = differential[1:] * differential[:-1] < 0
            turning_points_loc = np.array(np.where(turning_points))
        
            if turning_points_loc.any():
                assert len(turning_points_loc < 2), 'Found more than 1 turning point'
                tpk_est = peak_check_x[turning_points_loc].flatten()[0]
                max_response_est = peak_check_y[turning_points_loc].flatten()[0]
        
        return tpk_est, max_response_est
                                
    def plot(self, out_path, scale_type='standard', plot_residuals=False, hist_axes = False, fit_stats = True, convert_kelvin = False, plot_peak_uncertainty = False, plot_value_uncertainty =False, show_estimates=False, manual_axis_limits=False, title_extra = ''):
        #General function to sort out plot data and call the right plotting function
        #This is another hideous thing, I guess in the future I shouldn't try and write one function to do everything under the sun
    
        textdata = [self.R2, self.AIC, self.BIC] #Added to plot to show fit quality
        title = '{}: {}'.format(self.index, self.species_name) #Graph Title
        
        plt_x = self.temps
        plt_y = self.responses
        
        plt_x_curve = self.smooth_x
        plt_y_curve = self.smooth_y
        
        #Reformat data based on scale
        
        if convert_kelvin:
            plt_x = plt_x - self.Tref
            plt_x_curve = plt_x_curve - self.Tref
            temp_unit = 'C'
        else:
            temp_unit = 'K'
            
        if scale_type == 'log':
            if plt_x_curve.min() < 0:
                print('Warning: some x values are sub 0, x axis coerced to be positive')
                plt_x_curve -= plt_x_curve.min()
                plt_x -= plt_x.min()
            plt_x, plt_y, plt_x_curve, plt_y_curve = np.log(plt_x), np.log(plt_y), np.log(plt_x_curve), np.log(plt_y_curve)
                
        if scale_type == 'arrhenius':
            plt_y, plt_y_curve = np.log(plt_y), np.log(plt_y_curve)
            plt_x, plt_x_curve = 1 / plt_x, 1 / plt_x_curve
        
        #Get correct text data
        if self.model_name_short == 'LM':
            mid_text = 'Slope: {0[0]:.2f} \nIntercept: {0[1]:.2f}'.format([self.slope, self.intercept])
        else:
            mid_text = 'E:  {0[0]:.2f}\nB0:  {0[1]:.2f}'.format([self.final_E, self.final_B0])
        
        text_all = mid_text + '\nR2:  {0[0]:.2f}\nAIC: {0[1]:.2f}\nBIC: {0[2]:.2f}'.format(textdata)
        
        #Create output name
        sp_name = str(self.species_name).replace(' ', '_')     
        pattern = re.compile('[\W]+')
        path_adj = pattern.sub('', sp_name)  #remove non alphanumeric chars
        output_path = out_path + '/{}_{}'.format(self.index, path_adj) + title_extra + '.png'
        
        print('\tWriting: {}'.format(output_path))
        
        if hist_axes:
            self.plot2(plt_x, plt_y, plt_x_curve, plt_y_curve, text_all, title, scale_type, output_path, plot_residuals, fit_stats, temp_unit, plot_peak_uncertainty, plot_value_uncertainty, manual_axis_limits)
        else:
            self.plot1(plt_x, plt_y, plt_x_curve, plt_y_curve, text_all, title, scale_type, output_path, plot_residuals, fit_stats, temp_unit, plot_peak_uncertainty, plot_value_uncertainty, show_estimates, manual_axis_limits)
            
    def plot1(self, plt_x, plt_y, plt_x_curve, plt_y_curve, text_all, title, scale_type, output_path, plot_residuals, fit_stats, temp_unit, plot_peak_uncertainty, plot_value_uncertainty, show_estimates, manual_axis_limits):
        #Function to plot graph without histogram axes - less seaborn dependency.
      
        f = plt.figure()
        sns.set_style("ticks", {'axes.grid': True})
        
        ax = f.add_subplot(111)

        if manual_axis_limits:
            ax.set_ylim((0, max(plt_y) * 1.1))
            ax.set_xlim((min(plt_x) * 0.95 , max(plt_x) * 1.05))

        #Plot actual observed data as a scatterplot
        plt.plot(plt_x, plt_y, marker='o', linestyle='None', color='green', alpha=0.7, zorder=1)
        
        #Plot fitted curve
        plt.plot(plt_x_curve, plt_y_curve, marker='None', color='royalblue', linewidth=3, zorder=2)

        if show_estimates: #Show estimates for Tpk and THL if they exist. Try and catch is a bit hacky (like everything else) but works. 
            plt.plot(plt_x[1:], np.diff(plt_y), marker='None', color='purple', linewidth=3, zorder=2)
            try:
                plt.axvline(x=self.T_pk, color='red', zorder = 5)
                plt.axvline(x=self.T_H, color='lime', zorder = 5)
                plt.axvline(x=self.T_H_L, color='deepskyblue', zorder = 5)
            except:
                pass
            try:
                plt.axvline(x=self.final_T_pk, color='darkred', ls='--', zorder = 6)
                plt.axvline(x=self.final_T_H_L, color='darkblue', ls='--', zorder = 6)
                plt.axvline(x=self.final_T_H, color='darkgreen', ls='--', zorder = 6)
            except:
                pass

        plt.title(title, fontsize=14, fontweight='bold')
        
        if scale_type == 'log':
            plt.xlabel('log(Temperature) ({0})'.format(temp_unit))
            plt.ylabel('Log(' + self.trait + ')')
        elif scale_type == 'arrhenius':
            plt.xlabel('1 / Temperature ({0})'.format(temp_unit))
            plt.ylabel('Log(' + self.trait + ')')     
        else:
            plt.xlabel('Temperature ({0})'.format(temp_unit))
            plt.ylabel(self.trait)
            
        if fit_stats:
            plt.text(0.05, 0.85, text_all, ha='left', va='center', transform=ax.transAxes, color='darkslategrey')
         
        if self.uncertainty_estimated and plot_peak_uncertainty:
            tpk_est = getattr(self, 'tpk_est', None)
            response_est = getattr(self, 'max_response_est', None)

            tpk_uncertainty = getattr(self, 'tpk_uncertainty', 0) * 0.5
            response_uncertainty = getattr(self, 'response_uncertainty', 0) * 0.5
            
            vars = (tpk_est, response_est)
            if all(vars) and 'NA' not in vars:
                plt.errorbar(plt_x_curve[plt_y_curve.argmax()], self.max_response_est, xerr=(tpk_uncertainty / 2), yerr = (response_uncertainty / 2), 
                             ecolor='r', marker='x', zorder=3, mfc = 'r', ms = 3, elinewidth=1.5)        

        if plot_value_uncertainty:
            x_uncert = getattr(self, 'uncertainty_x', 0) * 0.5
            y_uncert = getattr(self, 'uncertainty_y', 0) * 0.5   
            max_visible_uncertainty = 100

            if scale_type == 'log':
                if isinstance(y_uncert, (np.ndarray, np.generic)):
                    y_uncert[y_uncert > 50] = max_visible_uncertainty
                    y_uncert[y_uncert < 0] = 0.1
                    y_uncert = np.log(y_uncert)

                if isinstance(x_uncert, (np.ndarray, np.generic)):
                    x_uncert[x_uncert > 50] = max_visible_uncertainty
                    x_uncert[x_uncert < 0] = 0.1
                    x_uncert = np.log(x_uncert)

            elif scale_type == 'arrhenius':
                if isinstance(y_uncert, (np.ndarray, np.generic)):
                    y_uncert[y_uncert > 50] = max_visible_uncertainty
                    y_uncert[y_uncert < 0] = 0.1
                    y_uncert = np.log(y_uncert)

                if isinstance(x_uncert, (np.ndarray, np.generic)):
                    x_uncert[x_uncert > 50] = max_visible_uncertainty
                    x_uncert = 1 / x_uncert
            else:
                if isinstance(y_uncert, (np.ndarray, np.generic)):
                    y_uncert[y_uncert > 50] = max_visible_uncertainty

                if isinstance(x_uncert, (np.ndarray, np.generic)):
                    x_uncert[x_uncert > 50] = max_visible_uncertainty

            plt.errorbar(plt_x, plt_y, xerr=x_uncert, yerr = y_uncert, marker='x', capsize=0, ls='none', fmt='', zorder=1, mfc = 'r', ecolor='r', elinewidth =0.5, ms = 1, alpha=0.6)   
                
        sns.despine() #Remove top and right border
        
        if plot_residuals: #create an inset plot with residuals
            if self.model_name_short == 'LM':
                residual_x = self.temps
                yvals = self.intercept + (residual_x * self.slope)
                residuals = yvals - self.responses
            elif self.model_name_short == 'BA':
                residual_x = self.temps
                yvals = np.exp(self.final_B0) * np.exp(-self.final_E / (residual_x * self.k)) 
                residuals = yvals - self.responses
            else:
                residuals = self.model.residual
                residual_x = plt_x
                
            ax2 = f.add_axes([.7, .65, .2, .2])
            ax2.plot(residual_x, residuals, color='royalblue', marker='o', linestyle='None', ms = 2)

            quad_mod = np.polyfit(residual_x, residuals, 2)
            xp = np.linspace(min(residual_x), max(residual_x), 100)
            poly = np.poly1d(quad_mod)

            ax2.plot(xp, poly(xp), marker='None', color='red', linewidth=1)

            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)

            ax2.grid(False)
            plt.title("Residuals")

        plt.savefig(output_path, bbox_inches='tight') 
        plt.close()
        
    def plot2(self, plt_x, plt_y, plt_x_curve, plt_y_curve, text_all, title, scale_type, output_path, plot_residuals, fit_stats, temp_unit, plot_peak_uncertainty, plot_value_uncertainty, manual_axis_limits):
        "Plots the graph with histogram axes, your mileage may vary..."
        
        #Scale works much better if defined manually
        if manual_axis_limits:
            _ylim = (0, plt_y.max() * 1.1)
        else:
            _ylim = None

        if scale_type == 'standard':
            max_x = int(np.ceil(max(plt_x_curve) / 10.0)) * 10
            min_x = int(np.floor(min(plt_x_curve) / 10.0)) * 10
        else:
            if scale_type == 'log':
                divisor = 100
            else:
                divisor = 20  
                addit_y = max(plt_y) / divisor 
                _ylim = (min(plt_y) - addit_y, max(plt_y) + 4 * addit_y) #Helps keep the data points away from the text
            addit_x = max(plt_x) / divisor
            max_x = max(plt_x) + addit_x
            min_x = min(plt_x) - addit_x
            
        df = pd.DataFrame({'x': plt_x, 'y': plt_y}, columns=["x", "y"])
        
        #plot the data and its distribution
        with sns.axes_style("white", {'axes.grid': True}):
            g = sns.jointplot(x="x", y="y", color='green', data=df,
                              stat_func=None,
                              ylim = _ylim,
                              xlim=(min_x, max_x), 
                              joint_kws=dict(alpha=0.5),
                              marginal_kws=dict(bins=20))
        
        #Add the fitted curve
        g.ax_joint.plot(plt_x_curve, plt_y_curve, marker='None', color='royalblue', linewidth=3)
        
        if fit_stats:
            g.ax_joint.text(0.05, 0.85, text_all, ha='left', va='center', transform=g.ax_joint.transAxes, color='darkslategrey')
        
        if scale_type == 'log':
            g.set_axis_labels('log(Temperature) ({0})'.format(temp_unit), 'Log(' + self.trait + ')')
        elif scale_type == 'arrhenius':
            g.set_axis_labels('1 / Temperature ({0})'.format(temp_unit), 'Log(' + self.trait + ')')
        else:
            g.set_axis_labels('Temperature ({0})'.format(temp_unit), self.trait)

        if plot_residuals: #create an inset plot with residuals
            if self.model_name_short == 'LM':
                residual_x = self.temps
                yvals = self.intercept + (residual_x * self.slope)
                residuals = yvals - self.responses
            elif self.model_name_short == 'BA':
                residual_x = self.temps
                yvals = np.exp(self.final_B0) * np.exp(-self.final_E / (residual_x * self.k)) 
                residuals = yvals - self.responses
            else:
                residuals = self.model.residual
                residual_x = plt_x
                
            sns.set_style("white", {'axes.grid': True})    
            ax2 = g.fig.add_axes([.845, .845, .13, .13])
            ax2.plot(residual_x, residuals, marker='None', color='royalblue', linewidth=1)
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)
            sns.despine()

        plt.savefig(output_path, bbox_inches='tight') 
        plt.close()
        
    def estimate_uncertainty(self, n=1000):
        assert self.model_name_short in ('SCH_F_TPK', 'SCH_F', 'SCH_TF', 'SCH_TF_TPK'), \
        'Covariance is only calculable for Schoolfield Sharpe models fitted via NLS'
        self.uncertainty_estimated = True

        tpk_values = [] #storage for tpk distribution
        max_response_values = [] #storage for max response distribution
        E_values = []

        params_names = self.model.var_names
        parameter_dict = self.model.params.valuesdict() #ordered dictionary of parameter mean values
        parameter_values = np.array([parameter_dict[i] for i in params_names])
        covar_mat = self.model.covar #matrix of covariances, order of parameter names is preserved

        #a few models seem to fail to have a matrix, not sure why - failed fits?
        #Some models also have infinite covariance which causes errors    
        if covar_mat is not None and np.all(np.isfinite(covar_mat)): 
            pars = np.random.multivariate_normal(parameter_values, covar_mat, n) #sample from bivariate distribution 1000 times

            for parameter_set in pars:
                valuesdict = {name: parameter_set[index] for index, name in enumerate(params_names)} #create a dictionary of parameter vals
                tpk, response_pk = self.differentiate_schoolfield(self.temps, valuesdict)

                if not isinstance(tpk, str) and not isinstance(response_pk, str):
                    tpk_values.append(tpk)
                    max_response_values.append(response_pk)
                    E_values.append(valuesdict['E'])

            return tpk_values, max_response_values, E_values 

        else:
            print('\n\t Warning: No covariance matrix\n')
            return np.nan, np.nan, np.nan
    
    def parameters_dict(self):
        "Returns a dictionary of the final parameters"
        
        final_B0 = getattr(self, 'final_B0', "NA")
        final_E = getattr(self, 'final_E', "NA")
        final_estimated_T_pk = getattr(self, 'tpk_est', "NA")
        final_max_response = getattr(self, 'max_response_est', "NA")
        final_E_D = getattr(self, 'final_E_D', "NA")
        final_E_D_L = getattr(self, 'final_E_D_L', "NA")
        final_T_H = getattr(self, 'final_T_H', "NA")
        final_T_H_L = getattr(self, 'final_T_H_L', "NA")
    
        param_dict = {
        "B0": final_B0,
        "E": final_E,
        "TPK": final_estimated_T_pk,
        "MaxResp": final_max_response,
        "ED": final_E_D,
        "EDL": final_E_D_L,
        "TH": final_T_H,
        "THL":  final_T_H_L}
        
        return param_dict
    
    def __lt__(self, other):
        "By defining less than we can directly compare models using max() to determine which is best"
        if self.AIC == other.AIC:
            return self.BIC > other.BIC
        else:
            return self.AIC > other.AIC

    def __gt__(self, other):
        "greater than"
        if self.AIC == other.AIC:
            return self.BIC < other.BIC
        else:
            return self.AIC < other.AIC

    def __le__(self, other):
        "less than or equal"
        if self.AIC == other.AIC:
            return self.BIC >= other.BIC
        else:
            return self.AIC >= other.AIC

    def __ge__(self, other):
        "less than or equal"
        if self.AIC == other.AIC:
            return self.BIC <= other.BIC
        else:
            return self.AIC <= other.AIC
            
    def __eq__(self, other):
        return self.AIC == other.AIC and self.BIC == other.BIC
