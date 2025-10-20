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

# create neural network surrogate model
from Data.Models.FcNet import FcNet


print("pytorch version", torch.__version__, "\n")

## parser arguments
parser = argparse.ArgumentParser(description='Load SgEncd Green Function for DownStream Tasks')
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
##############################################################################################
print('Network Architecture:', "\n", model)


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-1. generate matlab figdata for SgEncd Green function ...')

task1_meshdata = h5py.File('Data/task1_mesh_Gpts.mat')
SmpPts_fixedY1 = torch.from_numpy(np.float32(task1_meshdata['task1_Gpts_FixedY1'])) # data type: float32 
SmpPts_fixedY1 = SmpPts_fixedY1.to(device)
SmpPts_fixedY2 = torch.from_numpy(np.float32(task1_meshdata['task1_Gpts_FixedY2'])) # data type: float32 
SmpPts_fixedY2 = SmpPts_fixedY2.to(device)

SmpPts_x, SmpPts_y1, SmpPts_y2 = SmpPts_fixedY1[:,0:2], SmpPts_fixedY1[:,2:4], SmpPts_fixedY2[:,2:4]
SmpPts_z1 = torch.log( torch.norm((SmpPts_x - SmpPts_y1), dim=1, keepdim=True))
SmpPts_z2 = torch.log( torch.norm((SmpPts_x - SmpPts_y2), dim=1, keepdim=True))

SmpPts_fixedY1 = torch.cat((SmpPts_x, SmpPts_y1, SmpPts_z1), 1)
SmpPts_fixedY1_sym = torch.cat((SmpPts_y1, SmpPts_x, SmpPts_z1), 1)
SmpPts_fixedY2 = torch.cat((SmpPts_x, SmpPts_y2, SmpPts_z2), 1)
SmpPts_fixedY2_sym = torch.cat((SmpPts_y2, SmpPts_x, SmpPts_z2), 1)

# compute NN prediction
with torch.no_grad():        
    NN_G_fixedY1 = 0.5 * ( model(SmpPts_fixedY1) + model(SmpPts_fixedY1_sym) )
    NN_G_fixedY2 = 0.5 * ( model(SmpPts_fixedY2) + model(SmpPts_fixedY2_sym) )

io.savemat(os.path.join('Data', 'SgEncdGreen_FigMesh_Value.mat'), {'NN_G_fixedY1':np.array(np.float64(NN_G_fixedY1.cpu().detach())), 
                                                'NN_G_fixedY2':np.array(np.float64(NN_G_fixedY2.cpu().detach()))})

print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> task-2. compute solution using SgEncd Green function ...')

task2_quadraturedata = h5py.File('Data/task2_mesh_Qpts.mat')

SmpPts_Qpts = torch.from_numpy(np.float32(task2_quadraturedata['task2_mesh_Qpts'])) # data type: float32 
SmpPts_Qpts = SmpPts_Qpts.to(device)

SmpPts_x, SmpPts_y = SmpPts_Qpts[:,0:2], SmpPts_Qpts[:,2:4]
SmpPts_z = torch.log( torch.norm((SmpPts_x - SmpPts_y), dim=1, keepdim=True))

SmpPts_Qpts = torch.cat((SmpPts_x, SmpPts_y, SmpPts_z), 1)
SmpPts_Qpts_sym = torch.cat((SmpPts_y, SmpPts_x, SmpPts_z), 1)

batchsize = int(1e7)

# compute NN prediction
with torch.no_grad():        
    for i in range(int(len(SmpPts_Qpts)//batchsize)):
        
        NN_G_Qpts = 0.5 * ( model(SmpPts_Qpts[i*batchsize:(i+1)*batchsize,:]) + model(SmpPts_Qpts_sym[i*batchsize:(i+1)*batchsize,:]) )
   
        io.savemat(os.path.join('Data', 'SgEncdGreen_Qpts_Value_{0}.mat'.format(i)), {'NN_G_Qpts':np.array(NN_G_Qpts.cpu().detach())})

NN_G_Qpts = 0.5 * ( model(SmpPts_Qpts[(len(SmpPts_Qpts)//batchsize)*batchsize:,:]) + model(SmpPts_Qpts_sym[(len(SmpPts_Qpts)//batchsize)*batchsize:,:]) )

io.savemat(os.path.join('Data', 'SgEncdGreen_Qpts_Value_{0}.mat'.format(6)), {'NN_G_Qpts':np.array(NN_G_Qpts.cpu().detach())})

print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> task-3. solve eigenproblem of SgEncd Green function ...')

task3_FEM_grids = h5py.File('Data/task3_mesh_Gpts.mat')
SmpPts_diag = torch.from_numpy(np.float32(task3_FEM_grids['task3_diag_Qpts']))
SmpPts_Fem = torch.from_numpy(np.float32(task3_FEM_grids['task3_mesh_Gpts'])) 
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
    NN_G_diag_Qpts =    0.5 * ( model(SmpPts_diag) + model(SmpPts_diag_sym) )   
    NN_G_Fem = 0.5 * ( model(SmpPts_Fem) + model(SmpPts_Fem_sym) )   

io.savemat(os.path.join('Data', 'SgEncdGreen_Fem_Value.mat'), {'NN_G_diag_Qpts':np.array(NN_G_diag_Qpts.cpu().detach()), 'NN_G_Fem':np.array(NN_G_Fem.cpu().detach())})

