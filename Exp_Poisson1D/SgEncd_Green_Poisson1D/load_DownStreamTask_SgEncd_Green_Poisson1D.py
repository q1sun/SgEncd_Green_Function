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
parser.add_argument('-c', '--checkpoint', default='Data/Checkpoints/Poisson1D/simulation_0', type=str, metavar='PATH', help='path to save checkpoint')
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
print('===> task-1. generate matlab figdata for SgEncd Green function ...')

task1_meshdata = h5py.File('Data/task1_mesh_Gpts.mat')
SmpPts_FigMesh = torch.from_numpy(np.float64(task1_meshdata['task1_mesh_Gpts'])) # data type: float64 
SmpPts_FigMesh = SmpPts_FigMesh.to(device)

# generate mesh grids for the full xy domain
SmpPts_x, SmpPts_y = SmpPts_FigMesh[:,[0]], SmpPts_FigMesh[:,[1]]
SmpPts_z = torch.abs(SmpPts_x - SmpPts_y)
SmpPts_meshXY = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_meshXY_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# generate mesh grids for the full xz domain at a fixed y-pt
SmpPts_x, SmpPts_z = SmpPts_FigMesh[:,[0]], SmpPts_FigMesh[:,[1]]
SmpPts_y = 0.3 * torch.ones_like(SmpPts_y)
SmpPts_fixedY =  torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_fixedY_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# generate trace points at a fixed y-pt
SmpPts_x = torch.from_numpy(np.float64(task1_meshdata['mesh_x'])) # data type: float64 
SmpPts_x = SmpPts_x.to(device)
SmpPts_y = 0.3 * torch.ones_like(SmpPts_x)
SmpPts_z = torch.abs(SmpPts_x - SmpPts_y)
SmpPts_fixedY_trace =  torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_fixedY_trace_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# compute NN predicution
with torch.no_grad():        
    NN_G_meshXY = 0.5 * ( model(SmpPts_meshXY) + model(SmpPts_meshXY_sym) )
    NN_G_fixedY = 0.5 * ( model(SmpPts_fixedY) + model(SmpPts_fixedY_sym) )
    NN_G_fixedY_trace = 0.5 * ( model(SmpPts_fixedY_trace) + model(SmpPts_fixedY_trace_sym) )
    

io.savemat(os.path.join('Data', 'SgEncdGreen_FigMesh_Value.mat'), {'NN_G_meshXY':np.array(NN_G_meshXY.cpu().detach()), 
                                                'NN_G_fixedY':np.array(NN_G_fixedY.cpu().detach()), 'NN_G_fixedY_trace':np.array(NN_G_fixedY_trace.cpu().detach())})

print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-2. compute solution using deep Green function ...')

task2_quadraturedata = h5py.File('Data/task2_mesh_Qpts.mat')

SmpPts_Qpts = torch.from_numpy(np.float64(task2_quadraturedata['task2_mesh_Qpts'])) # data type: float64 
SmpPts_Qpts = SmpPts_Qpts.to(device)

SmpPts_x, SmpPts_y = SmpPts_Qpts[:,[0]], SmpPts_Qpts[:,[1]]
SmpPts_z = torch.abs(SmpPts_x - SmpPts_y)
SmpPts_Qpts = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_Qpts_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# compute NN prediction
with torch.no_grad():        
    NN_G_Qpts = 0.5 * ( model(SmpPts_Qpts) + model(SmpPts_Qpts_sym) )

io.savemat(os.path.join('Data', 'SgEncdGreen_Qpts_Value.mat'), {'NN_G_Qpts':np.array(NN_G_Qpts.cpu().detach())})

print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> task-3. solve eigenproblem of SgEncd Green function ...')

task3_FEM_grids = h5py.File('Data/task3_mesh_Gpts.mat')
SmpPts_Fem = torch.from_numpy(np.float64(task3_FEM_grids['task3_mesh_Gpts'])) # data type: float64 
SmpPts_Fem = SmpPts_Fem.to(device)

SmpPts_x, SmpPts_y = SmpPts_Fem[:,[0]], SmpPts_Fem[:,[1]]
SmpPts_z = torch.abs(SmpPts_x - SmpPts_y)
SmpPts_Fem = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_Fem_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

# compute NN predicution
with torch.no_grad():        
    NN_G_Fem = 0.5 * ( model(SmpPts_Fem) + model(SmpPts_Fem_sym) )   

io.savemat(os.path.join('Data', 'SgEncdGreen_Fem_Value.mat'), {'NN_G_Fem':np.array(NN_G_Fem.cpu().detach())})

print('*', '-' * 45, '*', "\n", "\n")

##############################################################################################




