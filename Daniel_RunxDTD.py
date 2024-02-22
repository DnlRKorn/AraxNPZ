#!/usr/bin/env python
# coding: utf-8

# In[22]:


import os, sys
import pickle
import joblib
import numpy as np
from tqdm import tqdm

path_list = os.getcwd().split('/')
index = path_list.index('KGML-xDTD')
model_path = '/'.join(path_list[:(index+1)] + ['model_evaluation','models'])

def load_drp_module(model_path: str):
    file_path = os.path.join(model_path,'kgml_xdtd','drp_module','model.pt')
    fitModel = joblib.load(file_path)
    return fitModel


# In[2]:


drp_module = load_drp_module(model_path)


# In[3]:


mondo_embs = np.load('mondo_embeddings.npz')["arr_0"]
no_embs = np.load('normalized_drug_embeddings.npz')['arr_0']


# In[4]:


#no_embs.shape[0]


# In[7]:


def predictForOneMONDO(mondo_idx):
    mondo_emb = mondo_embs[mondo_idx]
    X1 = np.tile(mondo_emb, (no_embs.shape[0], 1))
    X2 = no_embs
    X = np.hstack((X2,X1))
    res_temp = drp_module.predict_proba(X)
    return res_temp


# In[10]:


#castleman = mondo_embs[13864]
#adal = no_embs[2208]


# In[11]:


#X = np.hstack((adal,castle))
#X = X.reshape((1,1024))


# In[20]:


#drp_module.predict_proba(X)


# In[22]:


#drp_module.predict_proba(X)[0][1]


# In[8]:


#castle_res = predictForOneMONDO(13864)


# In[17]:


#castle_res[:,1].shape


# In[18]:


#mondo_embs.shape


# In[25]:


#z = np.zeros((mondo_embs.shape[0],no_embs.shape[0]))
#z[0] = castle_res[:,1]

run_num = int(sys.argv[1])
range_start = run_num * 1000
range_end = min(21660,(run_num+1) * 1000)
range_size = range_end - range_start

z = np.zeros((range_size,no_embs.shape[0]))

#for i in tqdm(range(21660)):
for j,i in tqdm(enumerate(range(range_start,range_end)),total=range_size):
    z[j] = predictForOneMONDO(i)[:,1]
#    if(i>3):break
np.savez(f'mondos_to_normalized_drugs_{run_num}.npz',z)

