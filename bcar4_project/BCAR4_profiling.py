#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 12:01:06 2022

@author: jacewebster

Creation of a BCAR4+ profile, based on TCGA microarry data,
using machine learning

"""

import pandas as pd
import numpy as np
from numpy import mean
from numpy import absolute
import shap
from sklearn import preprocessing

from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import LeaveOneOut

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, plot_confusion_matrix
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

###### DATA IMPORT AND BASIC FORMATTING
# Load starting TCGA datasets
microarray = pd.read_csv("../data/TCGA/microarray/AgilentG4502A_07_3_BRCA.tsv", sep="\t")
fusions = pd.read_csv("../data/TCGA/fusions/TCGA_fusions_annot.bedpe", sep="\t", header=None)

# Isolate BCAR4 fusions found in BRCA patients and get their patient IDs
bcar_fusions = fusions[(fusions[19]=="BRCA") & (fusions[6].str.contains("BCAR4"))].copy()

# Grab samples with no BCAR4 expression
exp = pd.read_csv("../data/TCGA/TCGA-BRCA.htseq_fpkm.roi.tsv", sep="\t")
texp = exp.transpose()
nobcar4 = texp[texp[1]==0].copy()

def slice(start=None, stop=None, step=1):
    return np.vectorize(lambda x: x[start:stop:step], otypes=[str])

nobcarcols = list(set(slice(start=0,stop=15)(nobcar4.index).flat))
nobcarstr = [str(i) for i in nobcarcols]
keepcols = []
for x in nobcarstr:
    if(x in microarray.columns):
        keepcols.append(x)
#keep_cols = list(set(slice(start=0, stop=15)(bcar_fusions[18]).flat))
#keep_cols.insert(0, 'sample')
# Keep only some samples, because their microarray data is not available...
# The first 4 are BCAR4 fusions. The rest are BCAR4-
# with 0 BCAR4 expression per RNA-Seq
usable_cols = keepcols.copy() # Samples with no bcar4 expresion

usable_cols = []
usable_cols.insert(0,'TCGA-BH-A0DG-01') # Adding fusion samples
usable_cols.insert(0,'TCGA-AO-A0JL-01')
#usable_cols.insert(0,'TCGA-C8-A1HL-01') # Keep out for validation testing
#usable_cols.insert(0,'TCGA-A1-A0SO-01') # Keep out for validation testing
usable_cols.insert(0, 'sample')
#usable_cols = ['sample', 
#               'TCGA-BH-A0DG-01',
#               'TCGA-AO-A0JL-01', 
#               'TCGA-C8-A1HL-01',
#               'TCGA-A1-A0SO-01'] + nobcarcols


fusion_microarray = microarray[usable_cols]

# Load Andy's microarray experiments
experiments = pd.read_csv("../data/Profiling/Microarray_results_All_genes.FromNicole.csv")
# Keep columns of interest, namely
# gene names, up/downregulation result, and microarray data
keep_cols = ['GeneName', 'TUH_B4siR271.status', 'TUH_B4siR268.status',
             'HME_ZCBOE.status', 'HME_LBOE.status', 'HME_LBMUT585.status',
             'HME_B4OE.status', 'TUH_B4siR271.logFC', 'TUH_B4siR268.logFC',
             'HME_ZCBOE.logFC', 'HME_LBOE.logFC', 'HME_LBMUT585.logFC',
             'HME_B4OE.logFC', 'HME_CONTROL_1', 'HME_CONTROL_2', 'HME_B4OE_1',
             'HME_B4OE_2', 'HME_LBMUT585_1', 'HME_LBMUT585_2', 'HME_LBOE_1',
             'HME_LBOE_2', 'HME_ZCBOE_1', 'HME_ZCBOE_2', 'TUH_CONTROL_1',
             'TUH_CONTROL_2', 'TUH_B4siR268_1', 'TUH_B4siR268_2', 'TUH_B4siR271_1',
             'TUH_B4siR271_2']
experiment_results = experiments[keep_cols]

# Merge datasets
df = pd.merge(experiment_results, fusion_microarray, left_on='GeneName', right_on='sample')
df = df.drop('sample',1) # Drop redundant column

# Important to note that BCAR4 is dropped here. 
# BCAR4 is not present on the TCGA's microarray.
# It may actually be best to not include BCAR4
# on the panel, since as a lncRNA, it is not
# available on most public microarray datasets.

# End of initial data import, save dataframe
df.to_csv('../data/Profiling/custom_dataframe.csv', index=False) # With all TCGA no BCAR4
df.to_csv('../data/Profiling/custom_dataframe.minimal.csv', index=False) # Only use TCGA fusions

###########

######## DATA WRANGLING and ADDITIONAL FORMATTING
# Load previously built dataframe
df = pd.read_csv('../data/Profiling/custom_dataframe.csv')
df = pd.read_csv('../data/Profiling/custom_dataframe.minimal.csv')

# Only keep genes if they were differentially expressed in
# >x experiments
df['changed'] = np.where(df['TUH_B4siR271.status']!='unchanged', 1, 0)
df['changed'] = np.where(df['TUH_B4siR268.status']!='unchanged',
                         df['changed'] + 1, df['changed'])
df['changed'] = np.where(df['HME_ZCBOE.status']!='unchanged',
                         df['changed'] + 1, df['changed'])
df['changed'] = np.where(df['HME_LBOE.status']!='unchanged',
                         df['changed'] + 1, df['changed'])
df['changed'] = np.where(df['HME_B4OE.status']!='unchanged',
                         df['changed'] + 1, df['changed'])
df['changed'] = np.where(df['HME_LBMUT585.status']!='unchanged',
                         df['changed'] + 1, df['changed'])

subdf = df[df['changed']>=4].copy() # Changed in 4 out of 6 experiments

# Standardize columns, as the cell data and the TCGA
# data have not been processed in the same way.
# Changes all expression columns to have a range of 0-1.
# Not sure if this is the best normalization method
# to use here, but it works..
#cols_to_normalize = subdf.columns[13:247]
cols_to_normalize = subdf.columns[13:163]
subdf[cols_to_normalize] = preprocessing.MinMaxScaler().fit_transform(subdf[cols_to_normalize])

# Transpose so each gene is a column, rows are samples
tdf = pd.DataFrame(subdf.transpose())
tdf.rename(columns=tdf.iloc[0], inplace=True)

# Drop unneeded rows
dropped_rows = tdf.index[0:13]
tdf_sampleonly = tdf.drop(dropped_rows)
tdf_sampleonly = tdf_sampleonly.drop('changed')


# Add BCAR4 status column for supervised algorithm
bcar4_status = [0,0,1,1,
                2, 2, 1, 1, # Don't use the LBMUT samples (2)
                1,1,1,1,
                0,0,0,0,
                1,1]
                #,1,1]  #Left off for validation testing
# Use when including all TCGA non BCAR4
zeros = [0] * (len(tdf_sampleonly) - len(bcar4_status))
status = bcar4_status + zeros
# else
status = bcar4_status
#status = bcar4_status
tdf_sampleonly['BCAR_Status'] = status

# Save this for future use
tdf_sampleonly.to_csv('../data/Profiling/normalized.csv')
tdf_sampleonly.to_csv('../data/Profiling/normalized.minimal.csv')

########



######## MODEL TRAINING
df = pd.read_csv('../data/Profiling/normalized.minimal.csv', index_col=0)

# Drop rows that won't be used for training
dropped_rows = ['HME_LBMUT585_1', 'HME_LBMUT585_2']

df = df.drop(dropped_rows)
df = df.dropna(axis=1)
df = df[keep_columns]

# Separate input values from target labels
X = df.iloc[:, :-1].values
Y = df.iloc[:, -1].values

# Just trying some random stuff here
# Best result from this is LR and KNC both giving score = 0.8 sd = 0.244
models = []
models.append(('LR', LogisticRegression(solver='liblinear', multi_class='ovr')))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNC', KNeighborsClassifier(algorithm='brute')))
# Evaluate each model
results = []
names = []
print("Name: Score (Standard Dev)")
for name, model in models:
    kfold = KFold(n_splits=10, random_state=1, shuffle=True)
    cv_results = cross_val_score(model, X, Y, cv=kfold)
    results.append(cv_results)
    names.append(name)
    print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))
    
## Other attempt using leave one out cross-validation
df = pd.read_csv('../data/Profiling/normalized.csv', index_col=0)
df = pd.read_csv('../data/Profiling/normalized.minimal.csv', index_col=0)
dropped_rows = ['HME_LBMUT585_1', 'HME_LBMUT585_2']

df = df.drop(dropped_rows)
df = df.dropna(axis=1)
#keep_cols_example = ["NUPR1", "CDC20", "GADD45A", "TNFRSF9", "BCAR_Status"]
model_data = df
#TESTING
keep_columns = ["BCAR_Status", "GDF15", "BMP5", "SYT4", "HAP1", "SEMA3C", 
                "GEM", "ZNF433", "IVL", "STK32A", "AIRE",
                "DNAH11", "RGS9BP", "RCOR2", "ERN2", "SLC37A3", "CGN", "MYCBPAP", "TREML1"]
model_data = df[keep_columns]
X = model_data.drop(columns="BCAR_Status")
Y = model_data["BCAR_Status"]
cv = LeaveOneOut()
model = RandomForestClassifier(random_state=42)
scores = cross_val_score(model, X, Y, scoring="neg_mean_absolute_error",
                         cv=cv, n_jobs=1)
mean(absolute(scores))

rf = RandomForestClassifier(random_state=42)

custom_back_order, custom_sorted_combo_backward = custom_sequential_backward(df, "BCAR_Status", rf)


def custom_sequential_backward(model_data, target, classifier):
    test_columns = list(model_data.columns)
    test_columns.remove(target)
    
    feature_combos = pd.DataFrame()

    backward_order = []
    
    while test_columns:
        cols = test_columns + [target]
        data = model_data[cols]
        X = data.drop(columns=target)
        Y = data[target]
        
        cv = LeaveOneOut()
        x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.30, random_state=13)

        classifier.fit(x_train, y_train)
        scores = cross_val_score(classifier, X, Y, scoring="neg_mean_absolute_error",
                                 cv=cv, n_jobs=1)
        score = mean(absolute(scores))
        test_instance ={'combo' : ', '.join(test_columns), 'score' : score, 'length' : len(test_columns)}
        feature_combos = feature_combos.append(test_instance, ignore_index=True)
        
        feature_importance = pd.Series(classifier.feature_importances_, index=x_train.columns).sort_values(ascending=False)
        remove = feature_importance.index[-1]

        test_columns.remove(remove)
        backward_order.extend([remove])

    sorted_combos = feature_combos.sort_values('score', ascending=False).reset_index()
    
    return backward_order, sorted_combos
        
    

##### GOING CRAZY HERE
df = pd.read_csv('../data/Profiling/normalized.csv', index_col=0)
df= pd.read_csv('../data/Profiling/normalized.minimal.csv', index_col=0)

# Drop rows that won't be used for training
dropped_rows = ['HME_LBMUT585_1', 'HME_LBMUT585_2']
df = df.drop(dropped_rows)
df = df.dropna(axis=1)

rf_classifier = RandomForestClassifier(random_state=42)

backward_order, sorted_combos_backward = sequential_backward(df, 'BCAR_Status', rf_classifier)
forward_order, sorted_combos_forward = sequential_forward(df, 'BCAR_Status', rf_classifier)


#Example
rf_classifier = RandomForestClassifier(random_state=42)
#The below gives an accuracy of 1 when used just on microarray data from Nicole
# and the fusion positive TCGA samples. Others give accuracy = 1 on that small group
model_columns = ["BCAR_Status", "NUPR1", "CDC20", "GADD45A", "STK32A", "NR3C2", "PDIA4",
                 "CLDN1", "LBR", "ERN2", "ZNF44"]
#This one (and others) gave accuracy of 1 on all BCAR4=0 TCGA + TCGA fusions + experiments
model_columns = ["BCAR_Status", "PDIA4", "LBR", "CDC20", "GADD45A", "STK32A", "TNFRSF9"]

# TESTING
model_columns = ["BCAR_Status", "GDF15", "BMP5", "SYT4", "HAP1", "SEMA3C", 
                "GEM", "ZNF433", "IVL", "STK32A", "AIRE",
                "DNAH11", "RGS9BP", "RCOR2", "ERN2", "SLC37A3", "CGN", "MYCBPAP", "TREML1"]
model_data = df[model_columns]
X = model_data.drop(columns="BCAR_Status")
Y = model_data["BCAR_Status"]
x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.20, random_state=13)
rf_classifier.fit(x_train, y_train)
y_pred = rf_classifier.predict(x_test)
accuracy_score(y_test, y_pred, normalize = True)
explainer = shap.TreeExplainer(rf_classifier)
shap_values = explainer.shap_values(x_test)
shap.summary_plot(shap_values, x_test)

colors = [(1, 1, 1), ((106/256, 235/256, 245/256)), ((22/256, 159/256, 169/256)), ((14/256, 95/256, 101/256))] # first color is black, last is red
cm = LinearSegmentedColormap.from_list("Custom", colors, N=30)
sns.set(font_scale=2.0)
plt.rcParams.update({'font.size': 18})

confusion(classifier=rf_classifier, x_test=x_test, y_test=y_test, y_pred=y_pred, model='Kickoff Return - Random Forest', cmap=cm)



######
def sequential_backward(model_data, target, classifier):
    """
    Accept a classification model and the full train/test data. Remove features one by one
    by determining the feature whose removal results in the highest accuracy score. Return
    a list of features in order of removal and a dataframe with all feature combinations
    and their accuracy scores, sorted by highest accuracy.
    """
    test_columns = list(model_data.columns)
    test_columns.remove(target)
    
    feature_combos = pd.DataFrame()

    backward_order = []

    while test_columns:
        cols = test_columns + [target]
        data = model_data[cols]

        X = data.drop(columns=target)
        Y = data[target]
        
        x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.30, random_state=13)

        classifier.fit(x_train, y_train)
        y_pred = classifier.predict(x_test)
        score = accuracy_score(y_test, y_pred, normalize = True)

        test_instance = {'combo' : ', '.join(test_columns), 'score' : score, 'length' : len(test_columns)}
        feature_combos = feature_combos.append(test_instance, ignore_index=True)

        feature_importance = pd.Series(classifier.feature_importances_, index=x_train.columns).sort_values(ascending=False)
        remove = feature_importance.index[-1]

        test_columns.remove(remove)
        backward_order.extend([remove])

    sorted_combos = feature_combos.sort_values('score', ascending=False).reset_index()
    
    return backward_order, sorted_combos

def sequential_forward(model_data, target, classifier):
    """
    Accept a classification model and the full train/test data. Find the single feature
    that provides the best classification on its own. Then add featurues one by one
    by determining the feature whose addition results in the highest accuracy score. Return
    a list of features in order of addition and a dataframe with all feature combinations
    and their accuracy scores, sorted by highest accuracy.
    """
    test_columns = list(model_data.columns)
    test_columns.remove(target)
    
    feature_combos = pd.DataFrame()

    forward_order = []
    i = 1

    while test_columns:
        for col in test_columns:
            cols = [target] + forward_order + [col]
            data = model_data[cols].copy()
            
            X = data.drop(columns=target)
            Y = data[target]
            #print(X.columns)
            x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.20, random_state=13)

            classifier.fit(x_train, y_train)
            y_pred = classifier.predict(x_test)

            score = accuracy_score(y_test, y_pred, normalize = True)

            test_instance = {'combo' : ', '.join(cols), 'score' : score, 'length' : i}
            feature_combos = feature_combos.append(test_instance, ignore_index=True)

        remove = feature_combos[feature_combos['length']==i].sort_values('score', ascending=False).reset_index().loc[0,'combo'].split(', ')[-1]

        test_columns.remove(remove)
        forward_order.extend([remove])
        i += 1

    sorted_combos = feature_combos.sort_values('score', ascending=False).reset_index()
    sorted_combos = sorted_combos.drop_duplicates(subset='length', keep='first')
    
    return forward_order, sorted_combos


def confusion(classifier, x_test, y_test, y_pred, model, cmap):
    """
    Plot a confusion matrix for a given classifier and its test data.
    """
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(18,8))
    fig.patch.set_alpha(0)

    norm_list = ['true', 'pred']
    titles = ['Normalized by Actual', 'Normalized by Prediction']
    labels = [1, 0]
    
    fig.suptitle(model + ' Confusion Matrices', fontsize=24, fontweight='bold')

    for ax, norm, title in zip(axes.flatten(), norm_list, titles):
        plot_confusion_matrix(classifier, 
                              x_test, 
                              y_test, 
                              ax=ax, 
                              cmap=cmap,
                              normalize=norm,
                              values_format='.1%',
                              labels=labels
                             )
        ax.set_title(title, pad=15, fontsize=20, fontweight='bold')
        ax.set_xlabel('Predicted Classification', labelpad=20, fontsize=20, fontweight='bold')
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=18)
        ax.set_yticklabels(labels, ha='right', fontsize=18)
        ax.set_ylabel('True Classification', labelpad=20, fontsize=20, fontweight='bold')
        ax.grid(None)
    
    plt.tight_layout(pad=2.2)  
    plt.show();






#########################################################
# Test


microarray = pd.read_csv("../data/TCGA/microarray/AgilentG4502A_07_3_BRCA.tsv", sep="\t")
microarray_t = pd.DataFrame(microarray.transpose())
microarray_t.rename(columns=microarray_t.iloc[0], inplace=True)
keep_columns = ["NUPR1", "CDC20", "GADD45A", "STK32A", "NR3C2", "PDIA4",
                 "CLDN1", "LBR", "ERN2", "ZNF44"]
keep_columns = ["GDF15", "BMP5", "SYT4", "HAP1", "SEMA3C", 
                "GEM", "ZNF433", "IVL", "STK32A", "AIRE",
                "DNAH11", "RGS9BP", "RCOR2", "ERN2", "SLC37A3", "CGN", "MYCBPAP", "TREML1"]
test_data = microarray_t[keep_columns].copy()
dropped_rows = test_data.index[0]
test_data = test_data.drop(dropped_rows)
#Normalize as before
test_data[keep_columns] = preprocessing.MinMaxScaler().fit_transform(test_data[keep_columns])
test_data.dropna(inplace=True)
test_pred = rf_classifier.predict(test_data)
#test_pred = model.predict(test_data)
test_pred = pd.DataFrame(data=test_pred,
                         index=test_data.index,
                         columns=["Prediction"])
test_pred.to_csv("../data/TCGA/microarray/pred.tsv")





