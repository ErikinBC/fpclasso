#########################################################
##### ------ SCRIPT TO PULL IN ICL DATASETS ------ ###### 

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--dir_output", default="home",dest="dir_output",
                  help="Directory to save output")
(options, args) = parser.parse_args()

# Assign
dir_output = options.dir_output
#dir_output="/home/erik/Documents/projects/FPC_Lasso/current/output/data"
print(dir_output)
# Change directory
import os
os.chdir(dir_output)
import pandas as pd
import numpy as np

from pmlb import fetch_data
from pmlb import classification_dataset_names, regression_dataset_names

n_class, n_reg = len(classification_dataset_names), len(regression_dataset_names)
# Create a data.frame to store
cn_reg = np.array(['dataset','n','p'])
df_reg = pd.DataFrame(index=range(n_reg), columns=range(len(cn_reg)))
df_reg.columns = cn_reg

cn_class = np.array(['dataset','n','p','class'])
df_class = pd.DataFrame(index=range(n_class), columns=range(len(cn_class)))
df_class.columns = cn_class

##################################################
##### ------- REGRESSION DATASETS --------- ######

jj = 0 
for ff in regression_dataset_names:
    jj += 1
    print(ff)
    tmp_X, tmp_y = fetch_data(ff,return_X_y=True)
    tmp_n, tmp_p = tmp_X.shape
    df_reg.loc[jj-1,'n'] = tmp_n
    df_reg.loc[jj-1,'p'] = tmp_p
    df_reg.loc[jj-1,'dataset'] = ff

# Only use datasets with < 10K observations and > 2 features
df_reg_sub = df_reg[np.array(df_reg.n < 1e4) & np.array(df_reg.p > 2)]

# Use only one of the "fri" datasets
type(df_reg_sub.dataset)
idx_fri = df_reg_sub.dataset.str.contains('fri_')
# Keep original plus one of the random ones
idx_fri_keep = np.concatenate((np.where(~idx_fri)[0],np.array([np.where(idx_fri)[0][0]])))
idx_fri_keep.sort()
# Subset again
df_reg_sub = df_reg_sub[~idx_fri]
df_reg_sub.reset_index(drop=True,inplace=True)

# Loop through and save datasets
jj = 0 
for ff in df_reg_sub.dataset:
    jj += 1
    print(ff)
    f1 = ff.split('_')
    f2 = 'reg_'+'_'.join(f1[1:])
    tmp_df = fetch_data(ff,return_X_y=False)
    tmp_df = tmp_df.reindex(columns=np.concatenate((np.array(['target']),np.setdiff1d(tmp_df.columns,'target'))))
    tmp_df.rename(columns={'target':'y'},inplace=True)
    tmp_df.to_csv(path_or_buf=dir_output + '/' + f2 + '.csv')

######################################################
##### ------- CLASSIFICATION DATASETS --------- ######

jj = 0 
for ff in classification_dataset_names:
    jj += 1
    print(ff)
    tmp_X, tmp_y = fetch_data(ff,return_X_y=True)
    tmp_class = len(np.unique(tmp_y))
    tmp_n, tmp_p = tmp_X.shape
    df_class.loc[jj-1,'n'] = tmp_n
    df_class.loc[jj-1,'p'] = tmp_p
    df_class.loc[jj-1,'dataset'] = ff
    df_class.loc[jj-1,'class'] = tmp_class

df_class_sub = df_class[np.array(df_class.n < 1e4) & np.array(df_class.p > 2) & np.array(df_class['class'] == 2)]
df_class_sub.reset_index(drop=True, inplace=True)

# Remove any duplicate names
tmp_names = pd.value_counts(np.array([x[0] for x in df_class_sub.dataset.str.split('\\_|\\-')]))
tmp_names = tmp_names[np.where(tmp_names > 1)[0]]
tmp_names = np.array(tmp_names.index)
# Get the index of the names not the list
tmp_names_index = df_class_sub.dataset.str.contains('|'.join(tmp_names))

df_class_sub = df_class_sub[~tmp_names_index]
df_class_sub.reset_index(drop=True, inplace=True)

classification_dataset_names_sub = np.array(df_class_sub.dataset)

jj = 0 
for ff in classification_dataset_names_sub:
    jj += 1
    print(ff)
    tmp_df = fetch_data(ff,return_X_y=False)
    tmp_df = tmp_df.reindex(columns=np.concatenate((np.array(['target']),np.setdiff1d(tmp_df.columns,'target'))))
    tmp_df.rename(columns={'target':'y'},inplace=True)
    tmp_df.to_csv(path_or_buf=dir_output + '/' + 'class_' + ff + '.csv')



