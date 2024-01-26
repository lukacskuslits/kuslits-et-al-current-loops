import torch
import torchvision
import torchvision.datasets as datasets
import numpy as np
import scipy
import scipy.io
from skimage.transform import resize


class DannDataLoader(torch.utils.data.Dataset):
    def __init__(self,length=8000):     
      self.SimInput = np.load('data/SimInput.npy')
      self.SimLabel = np.load('data/SimLabel.npy')
      self.RealInput = np.load('data/RealInput.npy')
      print(self.SimInput.shape)
      print(self.SimLabel.shape)
      print(self.RealInput.shape)

      self.length=length


    def __len__(self):
        return self.length


    def __getitem__(self, idx):
        ContainsNone=True
        while ContainsNone:
            NumpyIndex=np.random.randint(1,11)
            self.SimInput = np.load('data/SimInput_'+str(NumpyIndex)+'.npy')
            self.SimLabel = np.load('data/SimLabel_'+str(NumpyIndex)+'.npy')
            NotFound=True
            while NotFound:
                try:
                    SimIndex=np.random.randint(0,self.SimInput.shape[0])
                    SimInputSample = self.SimInput[SimIndex,:,:,:]
                    SimLabelSample = self.SimLabel[SimIndex,:,:,:]
                    NotFound=False
                except:
                    pass
            RealIndex=np.random.randint(0,self.RealInput.shape[0])
            RealInputSample = self.RealInput[RealIndex,:,:,:]
            ContainsNone=False
            if np.any(np.isnan(SimInputSample)):
                ContainsNone=True
            if np.any(np.isnan(SimLabelSample)):
                ContainsNone=True
            if np.any(np.isnan(RealInputSample)):
                ContainsNone=True


        return [ SimInputSample.astype(np.float32),  SimLabelSample.astype(np.float32), RealInputSample.astype(np.float32) ] 
   
   
