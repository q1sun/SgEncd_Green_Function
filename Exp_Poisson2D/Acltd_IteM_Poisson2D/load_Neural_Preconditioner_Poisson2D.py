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
parser.add_argument('-c', '--checkpoint', default='Data/Checkpoints/Poisson2D/simulation_0', type=str, metavar='PATH', help='path to save checkpoint')
# args, unknown = parser.parse_known_args()
args = parser.parse_args()
##############################################################################################


##############################################################################################
## ------------------------- ##
## network setting
## ------------------------- ##
dim_prob = 2
dim_in = 2 * dim_prob + 1
width = 40
depth = 3
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
model.double()
##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-1. generate neural preconditioners for preconditioning / hybrid item ...')

index_num = 4
grid_data = h5py.File('Data/mesh_Gpts_h{0}.mat'.format(index_num))

SmpPts_diag = torch.from_numpy(np.float64(grid_data['diag_Qpts']))
SmpPts_Fem = torch.from_numpy(np.float64(grid_data['mesh_Gpts'])) 
SmpPts_diag = SmpPts_diag.to(device)
SmpPts_Fem = SmpPts_Fem.to(device)

# generate sample points around the diagonal for numerical integration
SmpPts_diag_x, SmpPts_diag_y = SmpPts_diag[:,0:2], SmpPts_diag[:,2:4]
SmpPts_diag_z = torch.log( torch.norm((SmpPts_diag_x - SmpPts_diag_y), dim=1, keepdim=True))
SmpPts_diag = torch.cat((SmpPts_diag_x, SmpPts_diag_y, SmpPts_diag_z), 1)
SmpPts_diag_sym = torch.cat((SmpPts_diag_y, SmpPts_diag_x, SmpPts_diag_z), 1)

# generate sample points on the finite element mesh
SmpPts_x, SmpPts_y = SmpPts_Fem[:,0:2], SmpPts_Fem[:,2:4]
SmpPts_z = torch.log( torch.norm((SmpPts_x - SmpPts_y), dim=1, keepdim=True))
SmpPts_Fem = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_Fem_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# compute NN prediction
with torch.no_grad():   
    NN_G_diag_Qpts = 0.5 * ( model(SmpPts_diag) + model(SmpPts_diag_sym) )   
    NN_G_Fem = 0.5 * ( model(SmpPts_Fem) + model(SmpPts_Fem_sym) )   

io.savemat('Data/NN_Preconditioner_h{0}.mat'.format(index_num), {'NN_G_diag_Qpts':np.array(NN_G_diag_Qpts.cpu().detach()), 'NN_G_Fem':np.array(NN_G_Fem.cpu().detach())})

##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-2. generate neural preconditioners for multigrid ...')

index_num = 4
grid_data = h5py.File('App_MultiGrid/mesh_Gpts_h{0}.mat'.format(index_num))

SmpPts_diag = torch.from_numpy(np.float64(grid_data['diag_Qpts']))
SmpPts_Fem = torch.from_numpy(np.float64(grid_data['mesh_Gpts'])) 
SmpPts_diag = SmpPts_diag.to(device)
SmpPts_Fem = SmpPts_Fem.to(device)

# generate sample points around the diagonal for numerical integration
SmpPts_diag_x, SmpPts_diag_y = SmpPts_diag[:,0:2], SmpPts_diag[:,2:4]
SmpPts_diag_z = torch.log( torch.norm((SmpPts_diag_x - SmpPts_diag_y), dim=1, keepdim=True))
SmpPts_diag = torch.cat((SmpPts_diag_x, SmpPts_diag_y, SmpPts_diag_z), 1)
SmpPts_diag_sym = torch.cat((SmpPts_diag_y, SmpPts_diag_x, SmpPts_diag_z), 1)

# generate sample points on the finite element mesh
SmpPts_x, SmpPts_y = SmpPts_Fem[:,0:2], SmpPts_Fem[:,2:4]
SmpPts_z = torch.log( torch.norm((SmpPts_x - SmpPts_y), dim=1, keepdim=True))
SmpPts_Fem = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_Fem_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# compute NN prediction
with torch.no_grad():   
    NN_G_diag_Qpts = 0.5 * ( model(SmpPts_diag) + model(SmpPts_diag_sym) )   
    NN_G_Fem = 0.5 * ( model(SmpPts_Fem) + model(SmpPts_Fem_sym) )   

io.savemat('App_MultiGrid/NN_Preconditioner_h{0}.mat'.format(index_num), {'NN_G_diag_Qpts':np.array(NN_G_diag_Qpts.cpu().detach()), 'NN_G_Fem':np.array(NN_G_Fem.cpu().detach())})

##############################################################################################








