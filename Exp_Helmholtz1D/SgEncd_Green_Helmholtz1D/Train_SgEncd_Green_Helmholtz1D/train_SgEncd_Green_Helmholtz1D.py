##############################################################################################
import torch
import torch.nn as nn
import numpy as np
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import time
import datetime
import argparse 

from torch import optim, autograd
from matplotlib import pyplot as plt

# matlab to tensor and vice versa
from scipy.io import loadmat
import numpy as np
import scipy.io as io

# create training and testing datasets
from torch.utils.data import Dataset, DataLoader
from DataSets.Line1D import Sample_Points, Exact_Solution
from Utils import helper

# create neural network surrogate model
from Models.FcNet import FcNet

# load data from two datasets within the same loop
from itertools import cycle

print("pytorch version", torch.__version__, "\n")
torch.set_default_dtype(torch.float64)

## parser arguments
parser = argparse.ArgumentParser(description='SgEncd Green Function for 1D Helmholtz Equation using PINNs')
# checkpoints
parser.add_argument('-c', '--checkpoint', default='Checkpoints/Helmholtz1D/simulation_0', type=str, metavar='PATH', help='path to save checkpoint')
# datasets
parser.add_argument('-i', '--image', default='Images/Helmholtz1D/simulation_0', type=str, metavar='PATH', help='path to save figures')
args, unknown = parser.parse_known_args()
# args = parser.parse_args()
##############################################################################################

##############################################################################################
## hyperparameter configuration
## ------------------------- ##
# problem setting
dim_prob = 1
dim_in = 2 * dim_prob + 1
## ------------------------- ##
# dataset setting
num_reglr_pts_x = 500
num_reglr_pts_y = 500
num_snglr_pts_y = 500

# plot SgEncd green's function on test datasets
num_test_pts = 100

batch_num = 20
batchsize = num_reglr_pts_y * num_reglr_pts_y // batch_num
## ------------------------- ##
# network setting
width = 40
depth = 2
## ------------------------- ##
# optimization setting
beta_bndry = 400
beta_symtr = 400
beta_snglr = 400
num_epochs = 40000
milestones = 8000, 15000, 35000, 
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> preparing training and testing datasets ...')

# training dataset for sample points inside the domain
class TraindataReglr(Dataset):    
    def __init__(self, num_reglr_pts_x, num_reglr_pts_y): 
        
        # SmpPts_Reglr[ x-SmpPt, y-SmpPt, 0, 1 ], each x-SmpPt aligns with all y-SmpPts inside (0,1)
        self.SmpPts_Reglr = Sample_Points.SmpPts_Reglr_Line1D(num_reglr_pts_x, num_reglr_pts_y)                         
        self.k_Reglr_SmpPts = Exact_Solution.k_Exact_Line1D(self.SmpPts_Reglr[:,[0]])                          
    
    def __len__(self):
        return len(self.SmpPts_Reglr)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Reglr[idx]
        k_SmpPt = self.k_Reglr_SmpPts[idx]

        return [SmpPt, k_SmpPt]
    

class TraindataBndry(Dataset):
    def __init__(self, num_reglr_pts_y):
        self.SmpPts_Bndry = Sample_Points.SmpPts_Bndry_Line1D(num_reglr_pts_y)

    def __len__(self):
        return len(self.SmpPts_Bndry)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Bndry[idx]

        return SmpPt

class TraindataSnglr(Dataset):
    def __init__(self, num_snglr_pts):
        self.SmpPts_Snglr = Sample_Points.SmpPts_Snglr_Line1D(num_snglr_pts)

    def __len__(self):
        return len(self.SmpPts_Snglr)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Snglr[idx]

        return SmpPt
    
# testing dataset for equidistant sample points over the entire domain
class Testdata(Dataset):    
    def __init__(self, num_test_pts): 
        
        self.SmpPts_Test = Sample_Points.SmpPts_Test_Square2D(num_test_pts)   
              
    def __len__(self):
        return len(self.SmpPts_Test)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Test[idx]
 
        return SmpPt

# create training and testing datasets         
traindata_reglr = TraindataReglr(num_reglr_pts_x, num_reglr_pts_y)
traindata_bndry = TraindataBndry(num_reglr_pts_y)
traindata_snglr = TraindataSnglr(num_snglr_pts_y)
testdata = Testdata(num_test_pts)

# define dataloader 
dataloader_reglr = DataLoader(traindata_reglr, batch_size=batchsize, shuffle=True, num_workers=0)
dataloader_bndry = DataLoader(traindata_bndry, batch_size=batchsize, shuffle=True, num_workers=0)
dataloader_snglr = DataLoader(traindata_snglr, batch_size=batchsize, shuffle=True, num_workers=0)

dataloader_test = DataLoader(testdata, batch_size=num_test_pts*num_test_pts, shuffle=False, num_workers=0)

print('===> done!')
print('*', '-' * 45, '*')
##############################################################################################

##############################################################################################
# plot sample points during training and testing
if not os.path.isdir(args.image):
    helper.mkdir_p(args.image)
    
fig = plt.figure()
plt.scatter(traindata_reglr.SmpPts_Reglr[:,0], traindata_reglr.SmpPts_Reglr[:,1], c = 'red', label = 'interior points' )
plt.scatter(traindata_snglr.SmpPts_Snglr[:,0], traindata_snglr.SmpPts_Snglr[:,1], c = 'blue', label = 'interface points' )
plt.scatter(traindata_bndry.SmpPts_Bndry[:,0], traindata_bndry.SmpPts_Bndry[:,1], c = 'black', label = 'boundary points' )
plt.title('Sample Points during Training')
plt.legend(loc = 'lower right')
# plt.show()
fig.savefig(os.path.join(args.image,'TrainSmpPts.png'))

##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> creating training model ...')
print('*', '-' * 45, '*', "\n", "\n")

def train_epoch(epoch, model, optimizer, device):
    
    # set model to training mode
    model.train()

    training_loss_epoch, loss_reglr_epoch, loss_bndry_epoch, loss_symtry_epoch = 0, 0, 0, 0
    loss_snglr_epoch = 0

    for _, ( data_reglr, data_bndry, data_snglr ) in \
        enumerate(zip(dataloader_reglr, cycle(dataloader_bndry), cycle(dataloader_snglr)) ):
        
        # get mini-batch training data
        smppts_reglr, k_smppts = data_reglr
        smppts_bndry = data_bndry
        smppts_snglr = data_snglr

        # send training data to device            
        smppts_reglr = smppts_reglr.to(device)
        k_smppts = k_smppts.to(device)
        smppts_bndry = smppts_bndry.to(device)
        smppts_snglr = smppts_snglr.to(device)

        smppts_reglr_aug = torch.cat( (smppts_reglr, Exact_Solution.augmented_variable(smppts_reglr[:,[0]], smppts_reglr[:,[1]])), 1 )
        smppts_bndry_aug = torch.cat( (smppts_bndry, Exact_Solution.augmented_variable(smppts_bndry[:,[0]], smppts_bndry[:,[1]])), 1 )
        smppts_snglr_aug = torch.cat( (smppts_snglr, Exact_Solution.augmented_variable(smppts_snglr[:,[0]], smppts_snglr[:,[1]])), 1 )
        smppts_symtry_aug = torch.cat( (smppts_reglr[:,[1]], smppts_reglr[:,[0]], Exact_Solution.augmented_variable(smppts_reglr[:,[0]], smppts_reglr[:,[1]])), 1 )
       
        smppts_reglr_aug.requires_grad = True
        smppts_snglr_aug.requires_grad = True
        
        # forward pass to obtain NN prediction of correction function H(x,y)
        NN_reglr = model(smppts_reglr_aug)
        NN_bndry = model(smppts_bndry_aug)
        NN_snglr = model(smppts_snglr_aug)
        NN_symtry = model(smppts_symtry_aug)

        # zero parameter gradients and then compute NN prediction of {d^2/dy^2} H(x,y)
        model.zero_grad()
    
        grad_NN_reglr = torch.autograd.grad(outputs=NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0]
        grad_NN_snglr = torch.autograd.grad(outputs=NN_snglr, inputs=smppts_snglr_aug, grad_outputs=torch.ones_like(NN_snglr),retain_graph=True, create_graph=True, only_inputs=True)[0]
        Gx_NN_reglr, Gz_NN_reglr = grad_NN_reglr[:,[0]], grad_NN_reglr[:,[2]]
        Gz_NN_snglr = grad_NN_snglr[:,[2]]

        grad_Gx_NN_reglr = torch.autograd.grad(outputs=Gx_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gx_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0]
        grad_Gz_NN_reglr = torch.autograd.grad(outputs=Gz_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gz_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0]
        Gxx_NN_reglr, Gxz_NN_reglr = grad_Gx_NN_reglr[:,[0]], grad_Gx_NN_reglr[:,[2]]
        Gzz_NN_reglr = grad_Gz_NN_reglr[:,[2]]
        
        Residual_reglr = ( Exact_Solution.cx_Exact_Line1D(smppts_reglr[:,[0]]) * (Gx_NN_reglr + Gz_NN_reglr * (smppts_reglr[:,[0]] > smppts_reglr[:,[1]]) 
            - Gz_NN_reglr * (smppts_reglr[:,[0]] < smppts_reglr[:,[1]])) + Exact_Solution.c_Exact_Line1D(smppts_reglr[:,[0]]) * (Gxx_NN_reglr + Gzz_NN_reglr
            + 2 * (Gxz_NN_reglr * (smppts_reglr[:,[0]] > smppts_reglr[:,[1]]) - Gxz_NN_reglr * (smppts_reglr[:,[0]] < smppts_reglr[:,[1]]) ))
            + k_smppts ** 2 * NN_reglr )

        # construct mini-batch loss function and then perform backward pass
        loss_reglr = torch.mean( torch.pow( Residual_reglr, 2 ) )
        loss_snglr = torch.mean( torch.pow( 2 * Exact_Solution.c_Exact_Line1D(smppts_snglr[:,[0]]) * Gz_NN_snglr + torch.ones_like(Gz_NN_snglr), 2 ) )
        loss_bndry = torch.mean( torch.pow( NN_bndry, 2 ) )
        loss_symtry = torch.mean( torch.pow( NN_reglr - NN_symtry , 2) )
                                       
        loss_minibatch = loss_reglr + beta_bndry * loss_bndry + beta_snglr * loss_snglr + beta_symtr * loss_symtry

        # zero parameter gradients
        optimizer.zero_grad()
        # backpropagation
        loss_minibatch.backward()
        # parameter update
        optimizer.step()     

        # integrate loss over the entire training datset
        loss_reglr_epoch += loss_reglr.item() * smppts_reglr.size(0) / traindata_reglr.SmpPts_Reglr.shape[0]
        loss_bndry_epoch += loss_bndry.item() * smppts_bndry.size(0) / traindata_bndry.SmpPts_Bndry.shape[0]
        loss_snglr_epoch += loss_snglr.item() * smppts_snglr.size(0) / traindata_snglr.SmpPts_Snglr.shape[0]
        loss_symtry_epoch += loss_symtry.item() * smppts_reglr.size(0) / traindata_reglr.SmpPts_Reglr.shape[0]
        
        training_loss_epoch += loss_reglr_epoch + beta_bndry * loss_bndry_epoch + beta_snglr * loss_snglr_epoch + beta_symtr * loss_symtry_epoch

    return training_loss_epoch, loss_reglr_epoch, loss_bndry_epoch, loss_snglr_epoch, loss_symtry_epoch
##############################################################################################


##############################################################################################
print('*', '-' * 45, '*')
print('===> training neural network ...')

if not os.path.isdir(args.checkpoint):
    helper.mkdir_p(args.checkpoint)

# create model
model = FcNet.FcNet(dim_in, width, 1, depth)
model.Xavier_initi()
print('Network Architecture:', "\n", model)
print('Total number of trainable parameters = ', sum(p.numel() for p in model.parameters() if p.requires_grad))

# create optimizer and learning rate schedular
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, betas=(0.9, 0.999), eps=1e-08, weight_decay=0.01, amsgrad=False)
schedular = torch.optim.lr_scheduler.MultiStepLR(optimizer, milestones, gamma=0.1)

# load model to device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('DEVICE: {}'.format(device), "\n")
model = model.to(device)

# create log file
logger = helper.Logger(os.path.join(args.checkpoint, 'log.txt'), title='SgEncdGreen-Helmholtz-Line1D')
logger.set_names(['Learning Rate', 'TrainLoss', 'TrainL(reglr)', 'TrainL(bndry)', 'TrainL(snglr)', 'TrainL(symtry)'])
     
# train and test 
trainloss_curve, testloss_curve = [], []
trainloss_best = 1e10
since = time.time()
for epoch in range(num_epochs):
      
    print('Epoch {}/{}'.format(epoch, num_epochs-1), 'with LR = {:.1e}'.format(optimizer.param_groups[0]['lr']))  
        
    # execute training and testing
    trainloss_epoch, trainloss_reglr_epoch, trainloss_bndry_epoch, trainloss_snglr_epoch, trainloss_symtry_epoch = train_epoch(epoch, model, optimizer, device)
    
    # save current and best models to checkpoint
    is_best = trainloss_epoch < trainloss_best
    if is_best:
        print('==> Saving best model ...')
    trainloss_best = min(trainloss_epoch, trainloss_best)
    helper.save_checkpoint({'epoch': epoch + 1,
                            'state_dict': model.state_dict(),
                            'trainloss_reglr_epoch': trainloss_reglr_epoch,
                            'trainloss_bndry_epoch': trainloss_bndry_epoch,
                            'trainloss_snglr_epoch': trainloss_snglr_epoch,
                            'trainloss_symtry_epoch': trainloss_symtry_epoch,
                            'trainloss_epoch': trainloss_epoch,
                            'trainloss_best': trainloss_best,
                            'optimizer': optimizer.state_dict(),
                           }, is_best, checkpoint=args.checkpoint)   
    # save training process to log file
    logger.append([optimizer.param_groups[0]['lr'], trainloss_epoch, trainloss_reglr_epoch, 
                   trainloss_bndry_epoch, trainloss_snglr_epoch, trainloss_symtry_epoch])
    
    # adjust learning rate according to predefined schedule
    schedular.step()

    # print results
    trainloss_curve.append(trainloss_epoch)
    print('==> Full-Batch Training Loss = {:.4e}'.format(trainloss_epoch),"\n")
    
logger.close()
time_elapsed = time.time() - since

print('Done in {}'.format(str(datetime.timedelta(seconds=time_elapsed))), '!')
print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################

##############################################################################################
# plot learning curves
fig = plt.figure()
plt.plot(trainloss_curve, c = 'red', label = 'training loss' )
plt.title('Learning Curve during Training')
plt.legend(loc = 'upper right')
# plt.show()
fig.savefig(os.path.join(args.image,'TrainCurve.png'))
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> loading trained model for inference ...')

# load trained model
checkpoint = torch.load(os.path.join(args.checkpoint, 'model_best.pth.tar'))
model.load_state_dict(checkpoint['state_dict'])

# compute NN predicution of correction function H(x,y)
with torch.no_grad():    
    SmpPts_Test = testdata.SmpPts_Test
    SmpPts_Test = SmpPts_Test.to(device)
    NN_Test = model(torch.cat((SmpPts_Test,torch.abs(SmpPts_Test[:,0:1]-SmpPts_Test[:,1:2])),1)) 
    NN_Test = NN_Test.to('cpu')

# plot results
x = testdata.SmpPts_Test[:,0].reshape(num_test_pts,num_test_pts)
y = testdata.SmpPts_Test[:,1].reshape(num_test_pts,num_test_pts)

with torch.no_grad(): 

    plt.figure()
    G_NN = NN_Test.reshape(num_test_pts, num_test_pts)  
    plt.contourf(x, y, G_NN, 40, cmap = 'jet')
    plt.title('Deep Green Function on Test Dataset')
    bounds = torch.linspace(G_NN.min(), G_NN.max(), 5)
    plt.colorbar(ticks=bounds)
    # plt.show()  
    plt.savefig(os.path.join(args.image,'Green_NN_TestData.png'))
    plt.close(fig)

        
print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################




