# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:49:33 2023

@author: Minjun Son
"""

#%% Prepping the input and output data
import numpy as np
import pandas as pd
import random
import math
import matplotlib
from matplotlib import pyplot as plt
import tensorflow as tf
from sklearn.neural_network import MLPClassifier   # Load MLP from sklearn
from sklearn.metrics import accuracy_score   # Import accuracy score 
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Dropout
from sklearn.model_selection import KFold
from pathlib import Path 
from tensorflow.keras.layers import Input, Concatenate, BatchNormalization
# Print tensor flow version
print(tf.__version__)
# Preset Matplotlib figure sizes.
matplotlib.rcParams['figure.figsize'] = [9, 6]
#%%
###############################################################################
# Loading the raw data (insert filepaths)
df_ = pd.read_csv(r'/Users/minjun.son/Desktop/1-Analysis/GoI_info.csv')
df0 = pd.read_csv(r'/Users/minjun.son/Desktop/1-Analysis/completemat_20230803.csv')
###############################################################################

###############################################################################
# Get rid of the genes that does not have any ATAC score.
# Actaully this is not such a good idea. These do have information about GE as well
# ind=[];
# for aa in range(df0.shape[0]):
#     if np.isnan(df0.iloc[aa,2]) and np.isnan(df0.iloc[aa,65]):
#         ind = np.append(ind, aa)
# df1 = df0.drop(labels=ind, axis=0)
# df1 = df1.reset_index(drop=True)
# The original df0 has 182 columns

# So I'm just assigning df1 to be same as df0
df1 = pd.read_csv(r'/Users/minjun.son/Desktop/1-Analysis/completemat_20230803.csv')
###############################################################################

###############################################################################
# Correcting nan entries
for aa in range(df1.shape[0]):
    if np.isnan(df1.iloc[aa,2]):
        te=0*np.ones([1,21]); te[0,0]=1100;
        df1.iloc[aa,2:23] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,23]):
        te=0*np.ones([1,21]); te[0,0]=1100;
        df1.iloc[aa,23:44] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,44]):
        te=0*np.ones([1,21]); te[0,0]=1100;
        df1.iloc[aa,44:65] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,65]):
        te=0*np.ones([1,21]); te[0,0]=3300;
        df1.iloc[aa,65:86] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,86]):
        te=0*np.ones([1,21]); te[0,0]=3300;
        df1.iloc[aa,86:107] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,107]):
        te=0*np.ones([1,21]); te[0,0]=3300;
        df1.iloc[aa,107:128] = pd.DataFrame(te)
    if np.isnan(df1.iloc[aa,129]):
        df1.iloc[aa,129] = 3300
    if np.isnan(df1.iloc[aa,131]):
        df1.iloc[aa,131] = 3300
    if np.isnan(df1.iloc[aa,133]):
        df1.iloc[aa,133] = 3300
    # w=df1.iloc[aa,0]; i=w.find('.'); w=w[0:i]    
    # ind=[x for x, val in enumerate(np.array(df_.iloc[:,0])) if val == w]
    # if df_.iloc[ind[0],2] == '-':
    #     if df1.iloc[aa,2] != 0: df1.iloc[aa,2] = -df1.iloc[aa,2];
    #     if df1.iloc[aa,23] != 0: df1.iloc[aa,23] = -df1.iloc[aa,23];
    #     if df1.iloc[aa,44] != 0: df1.iloc[aa,44] = -df1.iloc[aa,44];
    #     if df1.iloc[aa,65] != 0: df1.iloc[aa,65] = -df1.iloc[aa,65];
    #     if df1.iloc[aa,86] != 0: df1.iloc[aa,86] = -df1.iloc[aa,86];
    #     if df1.iloc[aa,107] != 0: df1.iloc[aa,107] = -df1.iloc[aa,107];
###############################################################################

###############################################################################
# Getting rid of FF samples
ind=[];
for aa in range(df1.shape[0]):
    w = df1.iloc[aa,0]; i=w.find('.FF');
    if i != -1:
        ind = np.append(ind, aa)
df1 = df1.drop(labels=ind, axis=0)
df1 = df1.reset_index(drop=True)        
###############################################################################

###############################################################################          
# Spliting dataframe in sample, gene expression (output), and motif and NFkB data (input)
sam=df1.iloc[:,0]; ge=df1.iloc[:,1]; in_ = df1.iloc[:, list(range(2, 135)) + [137, 143]]
###############################################################################

del aa, df_, te, ind, w, i

# #%% Grouping the gene expressions based on the level.
# def bgg(totl, inc, minsize, hm):
#     if hm == 2:
#        haha=np.empty((0,2));
#        for aa in list(np.arange(minsize, totl-minsize+.01, inc)):
#            haha=np.append(haha, np.array([[aa, totl-aa]]), axis=0);
#     if hm == 3:
#        haha=np.empty((0,3));
#        for aa in list(np.arange(minsize, totl-2*minsize+.01, inc)):
#            for bb in list(np.arange(aa+minsize, totl-minsize+.01, inc)):
#                haha=np.append(haha, np.array([[aa, bb-aa, totl-bb]]), axis=0);       
#     return haha

# ### Then, create the train and test data
# ### Here define following:
# ###     The distance between groups. Number of groups. Min size for each group. The increment. The range
# dist=.2; nog=4; minsize=.4; inc=.1; ler=[-1, -.1]; rer=[2, 4];
# ac=np.empty((0,2*nog))
# for aa in list(np.arange(ler[0], ler[1]+.01, inc)):
#     for bb in list(np.arange(rer[0], rer[1]+.01, inc)):
#         totl=bb-aa-dist*(nog-1)
#         temp=bgg(totl, inc, minsize, nog-2)
#         for cc in list(range(temp.shape[0])):
#             temp1=[float('-inf'),aa]
#             a1=aa;
#             for dd in list(temp[cc]):
#                 temp1.extend([a1+dist, a1+dist+dd])
#                 a1 = a1+dist+dd;
#             temp1.extend([bb, float('inf')])
#             ac=np.append(ac, np.array([temp1]), axis=0)

# del aa, bb, a1, cc, dd, dist, nog, minsize, inc, ler, rer, totl, temp, temp1    
#%%
###############################################################################
# ge_grp=[[float('-inf'), -.6], [-.2, .2], [.6, 1.6], [2, 3], [4, float('inf')]]
# ge_grp=[[float('-inf'), -.6], [-.3, .55], [1.6, 2.6], [3.1, float('inf')]]

# ge_grp=[[-4, -.6], [-.2, .5], [1.5, 2.5], [3.5, 4.5], [5.5, 13]] # This gives around 65%
ge_grp=[[float('-inf'), -.6], [-.2, .5], [1.6, 3.3], [4.3, float('inf')]] # This one gives around 71%

print(ge_grp)

# The input summary
# 0:21 -> r1 pk1 (Dist, CPM, CPMold, 18 TFs)
# 21:42 -> r1 pk2 (Dist, CPM, CPMold, 18 TFs)
# 42:63 -> r1 pk3 (Dist, CPM, CPMold, 18 TFs)
# 63:84 -> r2 pk1 (Dist, CPM, CPMold, 18 TFs)
# 84:105 -> r2 pk2 (Dist, CPM, CPMold, 18 TFs)
# 105:126 -> r2 pk3 (Dist, CPM, CPMold, 18 TFs)
# 126:132 -> EZH2 count/location, HDAC3 count/loc, NELFE count/loc
# 132:180 -> NFkB trace features

rtt=list(range(0,180))
# rtt=list(range(0,21)) + list(range(63,66)) + list(range(66,114))
###############################################################################

###############################################################################
# Group samples based on the gene expression levels
ind=[]; ind_c=[]; gr=-1*np.ones([len(ge),1])
for aa in range(len(ge_grp)):
    ind0=[x for x, val in enumerate(ge) if (val > ge_grp[aa][0] and val < ge_grp[aa][1])]
    ind.extend([ind0])
    ind_c.extend(ind0)
    gr[ind0]=aa;
    print('Number of genes in Gr.' + str(aa) + ': ' + str(len(ind0)))
###############################################################################

###############################################################################
# Only the selected genes in all groups will be processed.
sam1=sam.iloc[ind_c]; sam1=sam1.reset_index(drop=True);
ge1=ge.iloc[ind_c]; ge1=ge1.reset_index(drop=True);
in_1=in_.iloc[ind_c,:]; in_1=in_1.reset_index(drop=True);
gr1=gr[ind_c,:];
###############################################################################
# Split the input data into two datasets
in_1_part1 = in_1.iloc[:, 0:132]  # Columns 3-135 are 133 columns
in_1_part2 = in_1.iloc[:, [133, 134]]  # Columns 137 and 143 are the next two columns

in_1_part1 = np.array(in_1_part1)
in_1_part2 = np.array(in_1_part2)

# Define kfold split and model
kfold = KFold(n_splits=5, shuffle=True, random_state=42)

# Create the multi-input model using the functional API
def create_functional_multi_input_model(input_dim1, input_dim2, out_dim):
    # Define the first input
    input1 = Input(shape=(input_dim1,))
    x1 = Dense(132, activation='relu')(input1)
    x1 = Dropout(0.15)(x1)
    x1 = Dense(66, activation='relu')(x1)
    x1 = Dropout(0.15)(x1)
    
    # Define the second input
    input2 = Input(shape=(input_dim2,))
    x2 = Dense(16, activation='relu')(input2)
    x2 = Dropout(0.1)(x2)
    x2 = Dense(8, activation='relu')(x2)
    x2 = Dropout(0.1)(x2)
    
    # Concatenate the outputs of the two inputs
    combined = Concatenate()([x1, x2])
    
    # Add further layers on top of the concatenated outputs
    x = Dense(134, activation='relu')(combined)
    x = Dropout(0.15)(x)
    output = Dense(out_dim, activation='softmax')(x)

    
    # Define the multi-input model
    model = Model(inputs=[input1, input2], outputs=output)
    
    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    model.compile(loss=loss_fn, optimizer='adam', metrics=['accuracy'])
    
    return model

# This part normalizes for convergence and avoids overfitting
take_abs = 1  # take absolute value of distance. (without directionality)
log_mode = 0  # log mode is on when 1
for aa in range(in_1_part1.shape[1]):
    if take_abs == 1:
        in_1_part1[:, aa] = abs(in_1_part1[:, aa])
    else:
        tmin = in_1_part1[:, aa].min()
        in_1_part1[:, aa] = in_1_part1[:, aa] - tmin
    if log_mode == 1:
        in_1_part1[:, aa] = np.log(in_1_part1[:, aa] + 1)
    tmax = in_1_part1[:, aa].max()
    if tmax == 0:
        tmax = 1
    in_1_part1[:, aa] = in_1_part1[:, aa] / tmax

for aa in range(in_1_part2.shape[1]):
    if take_abs == 1:
        in_1_part2[:, aa] = abs(in_1_part2[:, aa])
    else:
        tmin = in_1_part2[:, aa].min()
        in_1_part2[:, aa] = in_1_part2[:, aa] - tmin
    if log_mode == 1:
        in_1_part2[:, aa] = np.log(in_1_part2[:, aa] + 1)
    tmax = in_1_part2[:, aa].max()
    if tmax == 0:
        tmax = 1
    in_1_part2[:, aa] = in_1_part2[:, aa] / tmax
###############################################################################

# To store accuracy values
accuracy_values = []
aa=1

###############################################################################
# Mode1: To do kfold cross validation
for train_index, test_index in kfold.split(in_1_part1):
    X_train1, X_test1 = in_1_part1[train_index], in_1_part1[test_index]
    X_train2, X_test2 = in_1_part2[train_index], in_1_part2[test_index]
    y_train, y_test = gr1[train_index], gr1[test_index]

    model = create_functional_multi_input_model(len(in_1_part1[0]), len(in_1_part2[0]), len(ge_grp))
    
    tf.random.set_seed(42)
    model.fit([X_train1, X_train2], y_train, epochs=50, batch_size=8, verbose=0)

    loss, accuracy = model.evaluate([X_test1, X_test2], y_test, verbose=2)
    accuracy_values.append(accuracy)
            
    y_pred = model.predict([X_test1, X_test2], verbose=0)
    y_pred = tf.argmax(y_pred, axis=1)
    
    ############ To store the results for each training
    df2 = pd.DataFrame()
    temp=list(sam1[test_index]);
    df2['genes'] = temp
    df2['gr.true'] = np.vstack(gr1[test_index])
    df2['Pred'] = np.array(y_pred)
    filepath = Path(r'/Users/minjun.son/Desktop/1-Analysis/comb3_pred_' + str(aa) + '.csv')
    df2.to_csv(filepath, index=False)
    aa=aa+1
    ###########################################################################

# Plotting the accuracy values as a bar plot
plt.figure(figsize=(10, 6))
plt.bar(range(1, 6), accuracy_values, color='skyblue')
plt.xlabel('Fold')
plt.ylabel('Accuracy')
plt.title('K-Fold Cross-Validation Accuracy')
plt.xticks(range(1, 6))
plt.ylim(0, 1)
plt.show()
print(sum(accuracy_values)/5)