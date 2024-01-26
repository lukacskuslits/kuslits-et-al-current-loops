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

import numpy as np
import cv2
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

import numpy as np
import cv2

import segmentation_models_pytorch as smp
from segmentation_models_pytorch.encoders import get_preprocessing_fn

#training parameters
UseCuda=True


Net = smp.Unet(encoder_name="densenet121", encoder_weights="imagenet", in_channels=2, classes=5,  activation=None ).cuda(1)

Net.load_state_dict(torch.load("models/unet++_300.pth"))
Net.eval()


RealInput = np.load('data/SimInput_5.npy')
RealOut = np.load('data/SimLabel_5.npy')
print(RealInput.shape)
#training parameters
UseCuda=True


def mapcriterion(out,label):
    return (out-label).square().mean()


for i in range(RealInput.shape[0]):
     RealInputSample=RealInput[i,:,:,:]
     RealOut=torch.tensor(RealOut).cuda(1)
     
     RealInputSample.astype(np.float32)
     np.save("sim_out/input_"+str(i)+".npy",RealInputSample)
     RealIn=torch.tensor(RealInputSample).unsqueeze(0)
     RealIn = RealIn.cuda(1).float()
     alpha=0.0
     realeresponse = Net(RealIn)
     print( mapcriterion(realeresponse, RealOut[i,:,:,:]).cpu().detach().numpy() )

     np.save("sim_out/output_"+str(i)+".npy",realeresponse.cpu().detach().numpy())
     
     np.save("sim_out/label_"+str(i)+".npy",RealOut[i,:,:,:].cpu().detach().numpy())
     cv2.imwrite("sim_out/simres_"+str(i)+".png",realeresponse[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
   
