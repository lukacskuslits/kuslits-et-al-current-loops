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

#network creation
class DoubleConv(nn.Module):
    """Double convolution block...contains two conv layers with batchnorm and RELU"""
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.double_conv = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True)
        )

    def forward(self, x):
        return self.double_conv(x)


class Down(nn.Module):
    """Downscaling block...one downscale/maxpool and double conv"""
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.maxpool_conv = nn.Sequential(
            nn.MaxPool2d(2),
            DoubleConv(in_channels, out_channels)
        )

    def forward(self, x):
        return self.maxpool_conv(x)


class Up(nn.Module):
    """Upscaling bloc...upscaling by bilinaer upsambling, then double conv"""
    def __init__(self, in_channels, out_channels, bilinear=True):
        super().__init__()

        # if bilinear, use the normal convolutions to reduce the number of channels
        if bilinear:
            self.up = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
        else:
            self.up = nn.ConvTranspose2d(in_channels // 2, in_channels // 2, kernel_size=2, stride=2)

        self.conv = DoubleConv(in_channels, out_channels)

    def forward(self, x1, x2):
        x1 = self.up(x1)
        # input is CHW
        diffY = torch.tensor([x2.size()[2] - x1.size()[2]])
        diffX = torch.tensor([x2.size()[3] - x1.size()[3]])

        x1 = F.pad(x1, [diffX // 2, diffX - diffX // 2,
                        diffY // 2, diffY - diffY // 2])
        # if you have padding issues, see
        # https://github.com/HaiyongJiang/U-Net-Pytorch-Unstructured-Buggy/commit/0e854509c2cea854e247a9c615f175f76fbb2e3a
        # https://github.com/xiaopeng-liao/Pytorch-UNet/commit/8ebac70e633bac59fc22bb5195e513d5832fb3bd
        x = torch.cat([x2, x1], dim=1)
        return self.conv(x)


class OutConv(nn.Module):
    """Last output convolution of the network"""
    def __init__(self, in_channels, out_channels):
        super(OutConv, self).__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size=1)
    def forward(self, x):
        return self.conv(x)

class GradReverse(Function):
    @staticmethod
    def forward(ctx, x, lamda):
        ctx.lamda = lamda
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        output = (grad_output.neg() * ctx.lamda)
        return output, None



class UNet(nn.Module):
    """creation of the architecture from the previous blocks"""
    def __init__(self, n_channels, n_classes, bilinear=True):
        super(UNet, self).__init__()
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.bilinear = bilinear

        self.inc = DoubleConv(n_channels, 64)
        self.down1 = Down(64, 128)
        self.down2 = Down(128, 256)
        self.down3 = Down(256, 512)
        self.down4 = Down(512, 512)
        
        self.up1 = Up(1024, 256, bilinear)
        self.up2 = Up(512, 128, bilinear)
        self.up3 = Up(256, 64, bilinear)
        self.up4 = Up(128, 32, bilinear)
        
        #domain classifier
        self.domainconv1 = nn.Conv2d(32, 8, kernel_size=3, stride=2)
        self.domainconv2 = nn.Conv2d(8, 4, kernel_size=3, stride=2)
        self.domainfc1 = nn.Linear(4*10*21, 100) 
        self.domainfc2 = nn.Linear(100, 1) 
        
        
        self.outc = OutConv(32, n_classes)

    def forward(self, x, alpha):
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)

        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        
        x = self.up3(x, x2)
        x = self.up4(x, x1)
        
        reverse_feature = GradReverse.apply(x, alpha)
        
        d=self.domainconv1 (reverse_feature)
        d = F.leaky_relu(d)
        d=self.domainconv2 (d)
        d = F.leaky_relu(d)
        
        d= d.view(-1,4*10*21)
        d=self.domainfc1 (d)
        d = F.leaky_relu(d)
        d = self.domainfc2(d)
                
        logits = self.outc(x)
        return F.sigmoid(logits),  F.sigmoid(d)




RealInput = np.load('data/RealInput.npy')
print(RealInput.shape)
#training parameters
UseCuda=True

device = torch.device("cuda" if UseCuda else "cpu")

Net=UNet(2,5).to(device)

Net.load_state_dict(torch.load("models/unet_20.pth"))
Net.eval()
 
for i in range(RealInput.shape[0]):
     RealInputSample=RealInput[i,:,:,:]
     print(RealInputSample[0,0,0])
     RealInputSample.astype(np.float32)
     np.save("real_out/input"+str(i)+".npy",RealInputSample)
     RealIn=torch.tensor(RealInputSample).unsqueeze(0)
     RealIn = RealIn.to(device).float()
     alpha=0.0
     realeresponse, _ = Net(RealIn,alpha)
     np.save("real_out/output"+str(i)+".npy",realeresponse.cpu().detach().numpy())
     cv2.imwrite("real_out/realres_"+str(i)+".png",realeresponse[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
   
