import pandas as pd
import numpy as np
from Thermal_Models import read_database, clean_paired_data, gmean, normalise, weighted_std, weighted_amean

print('\n\tAggregating Summary...\n')
    
def mesophile_thermophile(entry):
    if entry > 323.15:
        return 'Thermophile'
    else:
        return 'Mesophile' 
    
agg_data = read_database('../Data/summaries/summary.csv')
agg_data = agg_data[(agg_data.Citation != "Mackey et al. (2013)")] # I really don't trust this data having looked at it (it came from BioTraits, not us) - TS
colnames = list(agg_data.columns.values)

#columns to keep
keep = ['Species', 'Trait', 'ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'RespiratoryPathway', 'Max_response', 'Est_Tpk', 'Rank', 'Points_Before_Peak', 'Tpk_std', 'Response_std', 'E_std', 'E']

#Drop any column not in the keep list
for colname in colnames:
    if colname not in keep:
        agg_data.drop(colname, 1, inplace=True)
       
#select best fits for growth rates only
groups = agg_data[(agg_data.Rank == 1) & (agg_data.Trait != 'Specific Growth Rate') & (agg_data.Points_Before_Peak > 2)]

#Rename appropriately
groups.rename(columns={'Est_Tpk': 'Est_Tpk',
                       'Tpk_std': 'Est_Tpk_std',
                       'Max_response': 'Est_Response',
                       'Response_std': 'Est_Response_std', 
                       'E' : 'E', 
                       'E_std': 'E_std'}, inplace=True)

#sort alphabetically by species name
groups.index.name = 'Species'
groups.sort_index()

#Add a temperature preference column
groups['TempPref'] =  groups['Est_Tpk'].apply(mesophile_thermophile)

#Fill blanks
groups = groups.fillna('NA')

#Save csv
groups.to_csv('../Data/summaries/non_aggregate_data_fluxes.csv')   
