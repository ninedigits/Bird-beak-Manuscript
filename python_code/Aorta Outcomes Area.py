#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from functools import reduce


# In[2]:


def bold_significance(val):
    if type(val)==str:
        val = float(val.replace('<', ''))
    fontweight = 'bold' if val<0.05 else 'normal'
    return 'font-weight: {}'.format(fontweight)


# In[3]:


# Load in data
bb_data = pd.read_csv('bbdata_area.csv')
bb_data.set_index('pat_id',inplace=True)
bb_data = bb_data[bb_data.index!=15]

# Calculate effective diameter from area
bb_data['effective_diameter'] = bb_data['area'].apply(lambda x: np.sqrt((4 * x)/np.pi))

# Calculate percent oversizing
bb_data['percent_oversizing'] = 100*((bb_data['graft'] / bb_data['effective_diameter'])-1)

# Curvature * diameter metric
bb_data['CD'] = bb_data['effective_diameter'] * bb_data['curve']

# Group by BBH
bb_data['group'] = bb_data['bbh'].apply(lambda x: 'BB' if x >= 5 else 'NBB')

results_data_all = bb_data.copy()
results_data_all['group'] = 'All'
results_data = pd.concat([bb_data, results_data_all])


results_data.rename(columns={
                                'bbh': 'BBH (mm)', 
                                'curve':'Curvature (mm-1)', 
                                'effective_diameter':'Diameter (mm)', 
                                'percent_oversizing':'Graft Oversizing (%)',
                                'bba': 'BBA (deg)',
                                'bbl': 'BBL (mm)',
                                'graft': 'Proximal Graft Diameter (mm)',
                                'area': 'Aortic Area (mm2)'
    
                             }, inplace=True)
results_data_output = results_data.groupby('group').agg([np.mean, np.std]).reset_index()
results_data_output


# In[28]:


bb_data


# In[29]:


results_data_output
one_sigfig = ['BBL (mm)', 'BBH (mm)','Diameter (mm)','CD']
three_sigfig = ['Curvature (mm-1)']
results_temp = []
cols = ['Proximal Graft Diameter (mm)','Aortic Area (mm2)','Curvature (mm-1)', 'Diameter (mm)', 'CD', 'Graft Oversizing (%)','BBL (mm)', 'BBH (mm)', 'BBA (deg)']
for col in cols:
    fstring = "{0:.3f} ± {1:.3f}" if col in three_sigfig else ("{0:.1f} ± {1:.1f}" if col in one_sigfig else "{0:.0f} ± {1:.0f}")
    mean_std = results_data_output[col]
    mean_std['mean ± std'] = mean_std.apply(lambda x: fstring.format(x['mean'], x['std']),axis=1)
    mean_std.drop(columns=['mean', 'std'],inplace=True)
    mean_std.index = ['All', 'BB', 'NBB']
    mean_std.columns = pd.MultiIndex.from_product([[col], mean_std.columns])
    results_temp.append(mean_std)
results_data_output = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True),results_temp)
display(results_data_output)    
for feat in cols:
    bbgroup = results_data[results_data['group'] == 'BB']
    nobbgroup = results_data[results_data['group'] == 'NBB']
    s, p = stats.ttest_ind(bbgroup[feat], nobbgroup[feat])
    results_data_output.loc['p-value', feat] = ['{:.3f}'.format(p)]
results_data_output.fillna('',inplace=True)
    
results_data_output = results_data_output.transpose()
results_data_output = results_data_output.reindex(sorted(results_data_output.index.values))
results_data_output.style.applymap(bold_significance, subset='p-value')
#writer = pd.ExcelWriter("/Users/maxfrohlich/Dropbox/Stanford-SJSU-Manuscript/Manuscript/Figures/bb_groups_demo.xlsx")
#results_data_output.to_excel(writer, 'sheet1')


# In[4]:


bb_data
stats.pearsonr(bb_data['bbh'], bb_data['bbl'])


# In[7]:


results_data_corr = results_data[results_data.group == 'All']
bba_out = []
bbl_out = []
bbh_out = []
for feat in ['Proximal Graft Diameter (mm)', 'Aortic Area (mm2)', 'Curvature (mm-1)', 'Diameter (mm)', 'Graft Oversizing (%)', 'CD']:
    current_feat = results_data_corr[feat]
    bba = results_data_corr['BBA (deg)']
    bbl = results_data_corr['BBL (mm)']
    bbh = results_data_corr['BBH (mm)']
    bba_r, bba_p = stats.pearsonr(current_feat, bba)
    bbl_r, bbl_p = stats.pearsonr(current_feat, bbl)
    bbh_r, bbh_p = stats.pearsonr(current_feat, bbh)
    bba_out.append(pd.DataFrame({'r-value': [bba_r], 'p-value':[bba_p]},index=[feat]))
    bbl_out.append(pd.DataFrame({'r-value': [bbl_r], 'p-value':[bbl_p]},index=[feat]))
    bbh_out.append(pd.DataFrame({'r-value': [bbh_r], 'p-value':[bbh_p]},index=[feat]))
bbl_corr = pd.concat(bbl_out)
bba_corr = pd.concat(bba_out)
bbh_corr = pd.concat(bbh_out)
bbal_corr = pd.concat({'BBL Correlation':bbl_corr.sort_index(), 
                       'BBA Correlation': bba_corr.sort_index(),
                       'BBH Correlation': bbh_corr.sort_index()})

#bbal_corr.to_csv('/Users/maxfrohlich/Dropbox/figure_1_raw/bba_bbl.csv')


bbal_corr.style.applymap(bold_significance, subset='p-value')


# In[ ]:




