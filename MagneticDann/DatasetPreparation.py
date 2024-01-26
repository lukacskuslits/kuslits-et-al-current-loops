import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim
import torch.utils.data
import torchvision.transforms as transforms
import torchvision.datasets as datasets
from skimage.transform import resize

import scipy
import scipy.io
import numpy as np
import cv2
"""## Data Preparation"""

#DATASET IMPORT
# -----------------------------------------------------------------------------------------------------------
#read data
def normalize(mat):
    mean=np.nanmean(mat)
    mat=np.nan_to_num(mat,mean)
    mat-=np.amin(mat)
    mat/=np.amax(mat)
    return mat
   
#read data from the mat files and convert them to numpy arrays 
for i in range(1,11):
    mat1 = scipy.io.loadmat('simulated_data_new_'+str(i)+'.mat')
    print(mat1.keys())
    print(len(mat1['sim_data'][0,0]))

    maps=np.transpose(mat1['sim_data']['maps'][0,0], (2,0,1))
    svs=np.transpose(mat1['sim_data']['svs'][0,0], (2,0,1))
    pos=np.transpose(mat1['sim_data']['pos'][0,0], (2,0,1))
    depths=np.transpose(mat1['sim_data']['depths'][0,0], (2,0,1))
    dtIs=np.transpose(mat1['sim_data']['dtIs'][0,0], (2,0,1))
    rads=np.transpose(mat1['sim_data']['rads'][0,0], (2,0,1))
    Is=np.transpose(mat1['sim_data']['Is'][0,0], (2,0,1))
    
    print(maps.shape)
    maps = resize(maps, (maps.shape[0], 96, 192))
    print(maps.shape)
    svs = resize(svs, (svs.shape[0], 96, 192))
    pos = resize(pos, (pos.shape[0], 96, 192))
    depths = resize(depths, (depths.shape[0], 96, 192))
    dtIs = resize(dtIs, (dtIs.shape[0], 96, 192))
    rads = resize(rads, (rads.shape[0], 96, 192))
    Is = resize(Is, (Is.shape[0], 96, 192))
    

    maps=normalize(maps)
    svs=normalize(svs)
    pos=normalize(pos)
    depths=normalize(depths)
    dtIs=normalize(dtIs)
    rads=normalize(rads)
    Is=normalize(Is)


    np.save('data/sim_maps_'+str(i)+'.npy',maps)
    np.save('data/sim_svs_'+str(i)+'.npy',svs)
    np.save('data/sim_pos_'+str(i)+'.npy',pos)
    np.save('data/sim_depths_'+str(i)+'.npy',depths)
    np.save('data/sim_dtIs_'+str(i)+'.npy',dtIs)
    np.save('data/sim_rads_'+str(i)+'.npy',rads)
    np.save('data/sim_Is_'+str(i)+'.npy',Is)
    
    SimLabel=np.stack((pos,depths,dtIs,rads,Is),1)
    SimInput=np.stack((maps,svs),1)

    np.save('data/SimInput_'+str(i)+'.npy',SimInput)
    np.save('data/SimLabel_'+str(i)+'.npy',SimLabel)

mat2 = scipy.io.loadmat('true_data_new.mat')
print(mat2.keys())
realmaps=np.transpose(mat2['true_data']['maps'][0,0], (2,0,1))
realsvs=np.transpose(mat2['true_data']['svs'][0,0], (2,0,1))

realmaps = resize(realmaps, (realmaps.shape[0], 96, 192))
realsvs = resize(realsvs, (realsvs.shape[0], 96, 192))
    
    
realmaps=normalize(realmaps)
realsvs=normalize(realsvs)



np.save('data/real_maps.npy',realmaps)
np.save('data/real_svs.npy',realsvs)


RealInput=np.stack((realmaps,realsvs),1)
    
np.save('data/RealInput.npy',RealInput)


#('labels', 'O'), ('maps', 'O'), ('maps_norm', 'O'), ('logmaps', 'O'), ('svs', 'O'), ('svs_norm', 'O'), ('pos', 'O'), #('rads', 'O'), ('Is', 'O'), ('depths', 'O'), ('dtIs', 'O')])
#mat2 = scipy.io.loadmat('output_pos_uNet_10_35.mat')['output_pos']
#mat3 = scipy.io.loadmat('BR_true.mat')['BR_true']
#print(mat2.shape)
#print(mat3.shape)




"""
mat1 = mat1.transpose(2,0,1)
mat2 = mat2.transpose(2,0,1)

#normalize between 0 and 1/ both input and output
mat1-=np.amin(mat1)
mat2-=np.amin(mat2)
mat1/=np.amax(mat1)
mat2/=np.amax(mat2)

#check intensity range
#print(np.amin(mat2))
#print(np.amax(mat2))

#reshape it to 1,1,45,90....batch, channel,width,height
mat1=mat1.reshape(8000,1,mat1.shape[1],mat1.shape[2])
mat2=mat2.reshape(8000,1,mat2.shape[1],mat2.shape[2])
print(mat1.shape)
print(mat2.shape)


# -----------------------------------------------------------------------------------------------------------
train_dataset = []
for i in range(1000):
    train_dataset.append([mat1[i,:,:], mat2[i,:,:]])
train_dataset = np.array(train_dataset)
validation_dataset = []
for i in range(1000,2000):
    validation_dataset.append([mat1[i,:,:], mat2[i,:,:]])
validation_dataset = np.array(validation_dataset)
target_dataset  = []
for i in range(2000,3000):
    target_dataset.append([mat1[i,:,:], mat2[i,:,:]])
target_dataset = np.array(train_dataset)
print(target_dataset.shape)
cv2.imwrite('example_0.png',(mat2[100,0,:,:]+1)*127.5)
cv2.imwrite('example_1000.png',(mat2[1100,0,:,:]+1)*127.5)
cv2.imwrite('example_2000.png',(mat2[2100,0,:,:]+1)*127.5)
cv2.imwrite('example_3000.png',(mat2[3100,0,:,:]+1)*127.5)
"""
