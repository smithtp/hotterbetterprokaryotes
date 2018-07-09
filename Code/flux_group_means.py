from Thermal_Models import read_database, clean_paired_data, gmean, normalise, weighted_std, weighted_amean
import numpy as np
import pandas as pd

print('\n\tCalculating EG values...\n')
    
    
agg_data = read_database('../Data/summaries/non_aggregate_data_fluxes_best.csv')

gb_cols = ['ConKingdom']

#calculate ES as the weighted geometric mean
def ES(g):
    pair_data = agg_data.ix[g.index]['E_std']
    return weighted_amean(g, pair_data, normalise_weights = True)

#calculate ES certainty as the weighted standard deviation 
def ES_min(g):
    pair_data = np.array(agg_data.ix[g.index]['E_std'])
    g = np.array(g)
    length = len(g)
    
    bootvals = []
    
    for i in range(1000):
        indices = np.random.choice(length, length) 
        g_boot = g[indices]
        pair_boot = pair_data[indices]
        bootvals.append(weighted_amean(g_boot, pair_boot, normalise_weights = True))
        
    return np.nanpercentile(bootvals, 2.5)

def ES_max(g):
    pair_data = np.array(agg_data.ix[g.index]['E_std'])
    g = np.array(g)
    length = len(g)
    
    bootvals = []
    
    for i in range(1000):
        indices = np.random.choice(length, length) 
        g_boot = g[indices]
        pair_boot = pair_data[indices]
        bootvals.append(weighted_amean(g_boot, pair_boot, normalise_weights = True))
        
    return np.nanpercentile(bootvals, 97.5)

for grouping in gb_cols:
    print(grouping)
    #group by each level
    group = agg_data.groupby(grouping)
    
    #aggregate
    group = group['E'].agg([len, ES, ES_max, ES_min])
    
    #rename the columns
    group.rename(columns={'len': grouping[0] + '_len',
                           'ES': grouping[0] + '_ES',
                           'ES_max': grouping[0] + '_ES_max',
                           'ES_min': grouping[0] + '_ES_min'}, inplace=True)
    
    #convert index to column so both aggdata and group have a matching column
    
    group.reset_index(level=0, inplace=True)
    if len(grouping) == 2:
        group.reset_index(inplace=True)
    
    #merge the group with the rest of the dataset
    agg_data = pd.merge(agg_data, group, how='outer', on=grouping)
    
agg_data.to_csv('../Data/summaries/group_means_non_aggregated_fluxes.csv')   
