##############################################################################################
import torch
import torch.nn as nn
import numpy as np
import os
import argparse 
import h5py

from torch import optim, autograd
from matplotlib import pyplot as plt

# matlab to tensor and vice versa
from scipy.io import loadmat
import numpy as np
import scipy.io as io
torch.set_default_dtype(torch.float64)

# create neural network surrogate model
from Data.Models.FcNet import FcNet

# load data from two datasets within the same loop
from itertools import cycle

print("pytorch version", torch.__version__, "\n")

## parser arguments
parser = argparse.ArgumentParser(description='Load SgEncd-Green Function for DownStream Tasks')
# checkpoints
parser.add_argument('-c', '--checkpoint', default='Data/Checkpoints/Helmholtz1D/simulation_0', type=str, metavar='PATH', help='path to save checkpoint')
# args, unknown = parser.parse_known_args()
args = parser.parse_args()
##############################################################################################


##############################################################################################
## ------------------------- ##
## network setting
## ------------------------- ##
dim_prob = 1
dim_in = 2 * dim_prob + 1
width = 40
depth = 2
## ------------------------- ##
## load model for downstream tasks
## ------------------------- ##
# create model
model = FcNet.FcNet(dim_in, width, 1, depth)
# load model to device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('DEVICE: {}'.format(device), "\n")
model = model.to(device)

# load trained model
checkpoint = torch.load(os.path.join(args.checkpoint, 'model_best.pth.tar'),map_location=torch.device('cpu'))
model.load_state_dict(checkpoint['state_dict'])
##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-1. generate neural preconditioners ...')

number_of_elements = 2 ** 12

quadraturedata = h5py.File('Data/mesh_Gpts_h{0}.mat'.format(number_of_elements))

SmpPts_Gpts = np.float64(quadraturedata['mesh_Gpts']) # data type: float64 
SmpPts_Gpts = torch.from_numpy(SmpPts_Gpts).t()

SmpPts_Gpts_Symtry = np.float64(quadraturedata['mesh_Gpts_symtry']) # data type: float64 
SmpPts_Gpts_Symtry = torch.from_numpy(SmpPts_Gpts_Symtry).t()

# compute NN predicution of Green's function
with torch.no_grad():    
    
    SmpPts_Gpts = SmpPts_Gpts.to(device)
    SmpPts_Gpts_Symtry = SmpPts_Gpts_Symtry.to(device)

    NN_Gpts = model(SmpPts_Gpts)
    NN_Gpts_Symtry = model(SmpPts_Gpts_Symtry)

    NN_Preconditioner = torch.reshape(0.5 * (NN_Gpts + NN_Gpts_Symtry), (number_of_elements - 1, number_of_elements - 1))    

io.savemat('Data/NN_Preconditioner_h{0}.mat'.format(number_of_elements), {'NN_Preconditioner':np.array(NN_Preconditioner.cpu().detach())})
##############################################################################################











