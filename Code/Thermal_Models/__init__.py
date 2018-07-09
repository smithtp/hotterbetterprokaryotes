#!/usr/bin/env python3

"""
This module provides a framework for and implementation of least squares and linear fitting of various 
thermal response models based on experimental data

Written in Python 3.5 Anaconda Distribution
"""

__title__ = 'Thermal Models'
__version__ = '0.0.1'
__author__ = 'jossthomas - jossthomas@maltbyhouse.co.uk'

#Anything imported into this file enters the name space of a script importing the module. Anything not imported is hidden. 

from .global_funcs import read_database, \
                          fit_models, \
                          split_datasets, \
                          rank_and_flatten, \
                          compile_models, \
                          bootstrap_model, \
                          estimate_uncertainty, \
                          clean_paired_data, \
                          gmean, \
                          normalise, \
                          weighted_std, \
                          weighted_amean, \
                          randomise_estimates, \
                          find_linear_arrhenius, \
                          plot_many_thermal, \
                          plot_weights_histograms

from .parameters import estimate_parameters

from .models import boltzmann_arrhenius, \
                    lm, \
                    schoolfield_two_factor, \
                    schoolfield_two_factor_tpk, \
                    schoolfield_full, \
                    schoolfield_full_tpk, \
                    boltzmann_arrhenius_two_weights_linear, \
                    boltzmann_arrhenius_two_weights_non_linear, \
                    boltzmann_arrhenius_two_weights_normalise_only, \
                    boltzmann_arrhenius_one_weight_linear, \
                    boltzmann_arrhenius_one_weight_non_linear
