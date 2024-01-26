# -*- coding: utf-8 -*-


"""# DANN training"""

import torch
import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim
import torch.nn.functional as F
import torch.utils.data
import torchvision.transforms as transforms
import torchvision.datasets as datasets
import danndataloader
from torch.autograd import Variable, Function
import segmentation_models_pytorch as smp
from segmentation_models_pytorch.encoders import get_preprocessing_fn

import numpy as np
import cv2
 
#train function
def train(Net, DataLoader, optimizer,  mapcritertion, domaincritertion, epoch,device):
    Net.train()
   
    start_steps = epoch * len(DataLoader)
    total_steps = epochs * len(DataLoader)
    for batch_idx, (SimIn, SimL, RealIn) in enumerate(DataLoader):
        SimIn = SimIn.to(device).float()
        SimL = SimL.to(device).float()
        RealIn = RealIn.to(device).float()
        
        optimizer.zero_grad()

        response = Net(SimIn)
        realeresponse = Net(RealIn)
        
     
        
        Classloss = mapcriterion(response, SimL)
        Loss=Classloss
        #loss= torch.mean(torch.abs(response-datay.float()))
        Loss.backward()
        #print loss in every iteration
        #print(loss.item())
        optimizer.step()
        if batch_idx%100==0:
                print("ClassLoss: "+str(Classloss.item())    )
    cv2.imwrite("dann/simres_"+str(epoch)+".png",response[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/simlabel_"+str(epoch)+".png",SimL[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/realres_"+str(epoch)+".png",realeresponse[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/simin_"+str(epoch)+".png",SimIn[0,0,:,:].cpu().detach().numpy()*255 )

#training parameters
epochs=500
UseCuda=True

batch_size = 32

train_loader = torch.utils.data.DataLoader(danndataloader.DannDataLoader(),batch_size=batch_size,
    drop_last=True,shuffle=True)

device = torch.device("cuda" if UseCuda else "cpu")

Net = smp.Unet(encoder_name="densenet121", encoder_weights="imagenet", in_channels=2, classes=5,  activation=None ).cuda()

def mapcriterion(out,label):
    return (out-label).square().mean()
    
domaincriterion = nn.BCELoss()
optimizer = torch.optim.Adam(Net.parameters())
for epoch in range(1, epochs):
        print("Epoch: "+str(epoch))
        train( Net, train_loader, optimizer,mapcriterion, domaincriterion, epoch, device)
        if (epoch%20==0):
            torch.save(Net.state_dict(), f"models/unet++_{epoch}.pth")



