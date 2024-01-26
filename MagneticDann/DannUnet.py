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
        self.domainfc1 = nn.Linear(4*10*22, 100) 
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
        d= torch.nn.functional.max_pool2d(d,2)
        
        
        d= d.view(-1,4*10*22)
        d=self.domainfc1 (d)
        d = F.leaky_relu(d)
        d = self.domainfc2(d)
                
        logits = self.outc(x)
        return logits,  F.sigmoid(d)


        
#train function
def train(Net, DataLoader, optimizer,  mapcritertion, domaincritertion, epoch,device):
    Net.train()
   
    start_steps = epoch * len(DataLoader)
    total_steps = epochs * len(DataLoader)
    for batch_idx, (SimIn, SimL, RealIn) in enumerate(DataLoader):
        SimIn = SimIn.to(device).float()
        SimL = SimL.to(device).float()
        RealIn = RealIn.to(device).float()
        
        p = float(batch_idx + start_steps) / total_steps
        alpha = 2. / (1. + np.exp(-10 * p)) - 1
        optimizer.zero_grad()
        response, sourcedomainpreds = Net(SimIn,alpha)
        realeresponse, targetdomainpreds = Net(RealIn,alpha)
        
     
        
        Classloss = mapcriterion(response, SimL)
        D_source = domaincritertion(sourcedomainpreds, torch.ones(sourcedomainpreds.shape[0], 1).to(device))
        D_target = domaincritertion(targetdomainpreds, torch.zeros(targetdomainpreds.shape[0], 1).to(device))
        DomainLoss = D_source+D_target
       
        Loss=Classloss+0.01*DomainLoss
        #Loss=DomainLoss
        #loss= torch.mean(torch.abs(response-datay.float()))
        Loss.backward()
        #print loss in every iteration
        #print(loss.item())
        optimizer.step()
        if batch_idx%100==0:
                print("ClassLoss: "+str(Classloss.item()) +" DomainLoss: "+str(DomainLoss.item())   )
    cv2.imwrite("dann/simres_"+str(epoch)+".png",response[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/simlabel_"+str(epoch)+".png",SimL[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/realres_"+str(epoch)+".png",realeresponse[0,:3,:,:].permute([1,2,0]).cpu().detach().numpy()*255 )
    cv2.imwrite("dann/simin_"+str(epoch)+".png",SimIn[0,0,:,:].cpu().detach().numpy()*255 )

#training parameters
epochs=500
UseCuda=True

batch_size = 8

train_loader = torch.utils.data.DataLoader(danndataloader.DannDataLoader(),batch_size=batch_size,
    drop_last=True,shuffle=True)

device = torch.device("cuda" if UseCuda else "cpu")

Net=UNet(2,5).to(device)

def mapcriterion(out,label):
    return (out-label).square().mean()
    
domaincriterion = nn.BCELoss()
optimizer = torch.optim.Adam(Net.parameters())
for epoch in range(1, epochs):
        print("Epoch: "+str(epoch))
        train( Net, train_loader, optimizer,mapcriterion, domaincriterion, epoch, device)
        if (epoch%20==0):
            torch.save(Net.state_dict(), f"models/unet_{epoch}.pth")



