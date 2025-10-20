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
from DataSets.Disc2D import Sample_Points, Exact_Solution
from Utils import helper

# create neural network surrogate model
from Models.FcNet import FcNet

# load data from two datasets within the same loop
from itertools import cycle

print("pytorch version", torch.__version__, "\n")
torch.set_default_dtype(torch.float64)

## parser arguments
parser = argparse.ArgumentParser(description='SgEncd Green Function for 2D Poisson Equation using PINNs')
# checkpoints
parser.add_argument('-c', '--checkpoint', default='Checkpoints/Poisson2D/simulation_0', type=str, metavar='PATH', help='path to save checkpoint')
# datasets
parser.add_argument('-i', '--image', default='Images/Poisson2D/simulation_0', type=str, metavar='PATH', help='path to save figures')
parser.add_argument('-r','--resume', default=False, type=bool, metavar='BOOL', help='whether to resume training from checkpoint')
args, unknown = parser.parse_known_args()
# args = parser.parse_args()
##############################################################################################

##############################################################################################
## hyperparameter configuration
## ------------------------- ##
# problem setting
dim_prob = 2
dim_aug = 1
dim_in = 2 * dim_prob + dim_aug
dim_out = 1
## ------------------------- ##
# dataset setting 
num_reglr_pts_x = 300
num_bndry_pts_x = 300
num_snglr_pts_x = 200

num_reglr_pts_y = 160

num_circles = 4
num_test_pts = 100 

batch_num = 8
batchsize_reglr = num_reglr_pts_y * num_reglr_pts_x // batch_num
batchsize_bndry = num_reglr_pts_y * num_bndry_pts_x // batch_num
batchsize_snglr = num_reglr_pts_y * num_circles // batch_num
## ------------------------- ##
# network setting
width = 40
depth = 3
## ------------------------- ##
# optimization setting
beta_bndry = 400 # penalty coefficient for Dirichlet boundary
beta_symtr = 400 # penalty coefficient for symmetry of Green's function
beta_snglr = 400
num_epochs = 30000
milestones = 10000, 20000,

##############################################################################################
# torch.autograd.set_detect_anomaly(True)

##############################################################################################
print('*', '-' * 45, '*')
print('===> preparing training and testing datasets ...')

# training dataset for sample points inside the domain
class TraindataReglr(Dataset):    
    def __init__(self, num_reglr_pts_x, num_reglr_pts_y): 
        
        # SmpPts_Train[ x-SmpPt, y-SmpPt ], each y-SmpPt aligns with all x-SmpPts inside B(0,1)
        self.SmpPts_Train = Sample_Points.SmpPts_Reglr_Disc2D(num_reglr_pts_x, num_reglr_pts_y)
        
    def __len__(self):
        return len(self.SmpPts_Train)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Train[idx]

        return SmpPt
    

class TraindataBndry(Dataset):
    def __init__(self, num_bndry_pts_x, num_reglr_pts_y):
        self.SmpPts_Bndry = Sample_Points.SmpPts_Bndry_Disc2D(num_bndry_pts_x, num_reglr_pts_y)

    def __len__(self):
        return len(self.SmpPts_Bndry)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Bndry[idx]

        return SmpPt

class TraindataSnglr(Dataset):
    def __init__(self, num_snglr_pts_x, num_reglr_pts_y, num_circles):
        self.SmpPts_Snglr = Sample_Points.SmpPts_Snglr_Disc2D(num_snglr_pts_x, num_reglr_pts_y, num_circles)

    def __len__(self):
        return len(self.SmpPts_Snglr)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Snglr[idx]

        return SmpPt
    
# testing dataset for equidistant sample points over the entire domain
class Testdata(Dataset):    
    def __init__(self, num_test_pts): 
        
        self.SmpPts_Test = Sample_Points.SmpPts_Test_Disc2D(num_test_pts)
        self.G_Exact_SmpPts = Exact_Solution.Green_Exact_Disc2D(self.SmpPts_Test[:,0:dim_prob], self.SmpPts_Test[:,dim_prob:2*dim_prob])     
              
    def __len__(self):
        return len(self.SmpPts_Test)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Test[idx]
        G_Exact_SmpPt = self.G_Exact_SmpPts[idx]
 
        return [SmpPt, G_Exact_SmpPt] 
    
    
class Plotdata(Dataset):    
    def __init__(self, num_test_pts): 
        
        self.SmpPts_Plot = Sample_Points.SmpPts_Plot_Disc2D(num_test_pts)
        self.G_Exact_SmpPts = Exact_Solution.Green_Exact_Disc2D(self.SmpPts_Plot[:,0:dim_prob], self.SmpPts_Plot[:,dim_prob:2*dim_prob])      
              
    def __len__(self):
        return len(self.SmpPts_Plot)
    
    def __getitem__(self, idx):
        SmpPt = self.SmpPts_Plot[idx]
        G_Exact_SmpPt = self.G_Exact_SmpPts[idx]
 
        return [SmpPt, G_Exact_SmpPt] 

# create training and testing datasets         
traindata_reglr = TraindataReglr(num_reglr_pts_x, num_reglr_pts_y)
traindata_bndry = TraindataBndry(num_bndry_pts_x, num_reglr_pts_y)
traindata_snglr = TraindataSnglr(num_snglr_pts_x, num_reglr_pts_y, num_circles)
testdata = Testdata(num_test_pts)   
plotdata = Plotdata(num_test_pts)

# define dataloader 
dataloader_reglr = DataLoader(traindata_reglr, batch_size=batchsize_reglr, shuffle=True, num_workers=0)
dataloader_bndry = DataLoader(traindata_bndry, batch_size=batchsize_bndry, shuffle=True, num_workers=0)
dataloader_snglr = DataLoader(traindata_snglr, batch_size=batchsize_snglr, shuffle=True, num_workers=0)

dataloader_test = DataLoader(testdata, batch_size=num_test_pts*num_test_pts, shuffle=False, num_workers=0)
dataloader_plot = DataLoader(plotdata, batch_size=num_test_pts*num_test_pts, shuffle=False, num_workers=0)

print('===> done!')
print('*', '-' * 45, '*')
##############################################################################################

##############################################################################################
# plot sample points during training and testing
if not os.path.isdir(args.image):
    helper.mkdir_p(args.image)
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> creating training model ...')
print('*', '-' * 45, '*', "\n", "\n")

def train_epoch(epoch, model, optimizer, device):
    
    # set model to training mode
    model.train()

    loss_epoch, loss_reglr_epoch, loss_bndry_epoch, loss_snglr_epoch, loss_symtry_epoch = 0, 0, 0, 0, 0

    for _, ( data_reglr, data_bndry, data_snglr ) in \
        enumerate(zip(dataloader_reglr, cycle(dataloader_bndry), cycle(dataloader_snglr)) ):
        
        # get mini-batch training data
        smppts_reglr = data_reglr.to(device)
        smppts_bndry = data_bndry.to(device)
        smppts_snglr = data_snglr.to(device)
        smppts_snglr = torch.reshape(smppts_snglr, (-1, 4))

        smppts_reglr_aug = torch.cat( (smppts_reglr, Exact_Solution.augmented_func(smppts_reglr)), dim=1 )
        smppts_bndry_aug = torch.cat( (smppts_bndry, Exact_Solution.augmented_func(smppts_bndry)), dim=1 )
        smppts_snglr_aug = torch.cat( (smppts_snglr, Exact_Solution.augmented_func(smppts_snglr)), dim=1 )
        smppts_symtry_aug = torch.cat( (smppts_reglr[:, dim_prob:2*dim_prob], smppts_reglr[:, 0:dim_prob], Exact_Solution.augmented_func(smppts_reglr) ), dim=1 )      
        
        smppts_reglr_aug.requires_grad = True
        smppts_snglr_aug.requires_grad = True
        
        NN_reglr = model(smppts_reglr_aug)
        NN_bndry = model(smppts_bndry_aug)
        NN_snglr = model(smppts_snglr_aug)
        NN_symtry = model(smppts_symtry_aug)

        # zero parameter gradients and then compute NN prediction of the derivatives 
        model.zero_grad()
        grad_NN_reglr = torch.autograd.grad(outputs=NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0]
        Gx1_NN_reglr, Gx2_NN_reglr, Gz_NN_reglr = grad_NN_reglr[:,[0]], grad_NN_reglr[:,[1]], grad_NN_reglr[:,[-1]]
        
        Gx1x1_NN_reglr = torch.autograd.grad(outputs=Gx1_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gx1_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0][:,[0]]
        Gx2x2_NN_reglr = torch.autograd.grad(outputs=Gx2_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gx2_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0][:,[1]]
        Gx1z_NN_reglr = torch.autograd.grad(outputs=Gx1_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gx1_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0][:,[-1]]
        Gx2z_NN_reglr = torch.autograd.grad(outputs=Gx2_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gx2_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0][:,[-1]]
        Gzz_NN_reglr = torch.autograd.grad(outputs=Gz_NN_reglr, inputs=smppts_reglr_aug, grad_outputs=torch.ones_like(Gz_NN_reglr), retain_graph=True, create_graph=True, only_inputs=True)[0][:,[-1]]
        
        dist_xy_reglr = torch.norm(smppts_reglr[:,0:2] - smppts_reglr[:,2:4], dim=1, keepdim=True)
        normVec_reglr = (smppts_reglr[:,0:2] - smppts_reglr[:,2:4]) / dist_xy_reglr ** 2
        
        Residual_reglr = - ( Gx1x1_NN_reglr + Gx2x2_NN_reglr + 2 * Gx1z_NN_reglr * normVec_reglr[:,[0]]
                      + 2 * Gx2z_NN_reglr * normVec_reglr[:,[1]] + Gzz_NN_reglr / dist_xy_reglr ** 2)
        
        grad_NN_snglr = torch.autograd.grad(outputs=NN_snglr, inputs=smppts_snglr_aug, grad_outputs=torch.ones_like(NN_snglr), retain_graph=True, create_graph=True, only_inputs=True)[0]
        Gx1_NN_snglr, Gx2_NN_snglr, Gz_NN_snglr = grad_NN_snglr[:,[0]], grad_NN_snglr[:,[1]], grad_NN_snglr[:,[-1]]

        Gx1_NN_snglr = torch.reshape(Gx1_NN_snglr, (-1, num_snglr_pts_x))
        Gx2_NN_snglr = torch.reshape(Gx2_NN_snglr, (-1, num_snglr_pts_x))
        diff_x1_y1 = torch.reshape(smppts_snglr[:,[0]] - smppts_snglr[:,[2]], (-1, num_snglr_pts_x))
        diff_x2_y2 = torch.reshape(smppts_snglr[:,[1]] - smppts_snglr[:,[3]], (-1, num_snglr_pts_x))
        Gz_NN_snglr = torch.reshape(Gz_NN_snglr, (-1, num_snglr_pts_x))

        Integral_Snglr = - 2 * np.pi * torch.mean(Gx1_NN_snglr * diff_x1_y1 + Gx2_NN_snglr * diff_x2_y2 + Gz_NN_snglr, dim=1)

        # Integral_Snglr = - 2 * np.pi * torch.mean(torch.reshape((Gx1_NN_snglr * (smppts_snglr[:,[0]] - smppts_snglr[:,[2]]) + Gx2_NN_snglr * (smppts_snglr[:,[1]] - smppts_snglr[:,[3]]) + Gz_NN_snglr), 
        #                                         (-1, num_snglr_pts_x) ), dim=1)
        
        # construct mini-batch loss function and then perform backward pass
        loss_reglr = torch.mean( torch.pow( torch.squeeze( Residual_reglr + torch.zeros_like(Residual_reglr)), 2 ) )
        loss_snglr = torch.mean( torch.pow( torch.squeeze( Integral_Snglr - 1 ), 2 ) )
        loss_bndry = torch.mean( torch.pow( torch.squeeze( NN_bndry ), 2 ) )  
        loss_symtry = torch.mean( torch.pow( torch.squeeze( NN_reglr - NN_symtry ), 2) )
                                       
        loss_minibatch = loss_reglr + beta_bndry * loss_bndry + beta_snglr * loss_snglr + \
                            beta_symtr * loss_symtry 

        # zero parameter gradients
        optimizer.zero_grad()
        # backpropagation
        loss_minibatch.backward()
        # parameter update
        optimizer.step()     

        # integrate loss over the entire training datset
        loss_reglr_epoch += loss_reglr.item() * smppts_reglr.size(0) / traindata_reglr.SmpPts_Train.shape[0]
        loss_bndry_epoch += loss_bndry.item() * smppts_bndry.size(0) / traindata_bndry.SmpPts_Bndry.shape[0]
        loss_snglr_epoch += loss_snglr.item() * Integral_Snglr.size(0) / traindata_snglr.SmpPts_Snglr.shape[0]
        loss_symtry_epoch += loss_symtry.item() * smppts_reglr.size(0) / traindata_reglr.SmpPts_Train.shape[0]

        loss_epoch += loss_reglr_epoch + beta_bndry * loss_bndry_epoch + beta_snglr * loss_snglr_epoch + \
                        beta_symtr * loss_symtry_epoch

    return loss_epoch, loss_reglr_epoch, loss_bndry_epoch, loss_snglr_epoch, loss_symtry_epoch
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> creating testing model ...')
print('*', '-' * 45, '*', "\n", "\n")

def test_epoch(epoch, model, optimizer, device):
    
    # set model to testing mode
    model.eval()

    testloss_epoch = 0

    for smppts_test, G_exact_smppts in dataloader_test:
        
        # send inputs, outputs to device
        smppts_test = smppts_test.to(device)
        G_exact_smppts = G_exact_smppts.to(device)
        
        # forward pass and then compute loss function for approximating G(x,y)
        NN_smppts = model(torch.cat((smppts_test, Exact_Solution.augmented_func(smppts_test)), 1)) 
        loss_epoch = torch.mean( torch.pow( torch.squeeze(NN_smppts) - G_exact_smppts, 2 ) )                 
                
        # integrate loss      
        testloss_epoch += loss_epoch.item()          
    
    
    return testloss_epoch
##############################################################################################



##############################################################################################
print('*', '-' * 45, '*')
print('===> training neural network ...')

if not os.path.isdir(args.checkpoint):
    helper.mkdir_p(args.checkpoint)

# create model
model = FcNet.FcNet(dim_in, width, dim_out, depth)

print('Network Architecture:', "\n", model)
print('Total number of trainable parameters = ', sum(p.numel() for p in model.parameters() if p.requires_grad))

# create optimizer and learning rate schedular
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, betas=(0.9, 0.999), eps=1e-08, weight_decay=0.01, amsgrad=False)
schedular = torch.optim.lr_scheduler.MultiStepLR(optimizer, milestones, gamma=0.1)

# load model to device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('DEVICE: {}'.format(device), "\n")
# model = nn.DataParallel(model)
model = model.to(device)

# create log file
logger = helper.Logger(os.path.join(args.checkpoint, 'log.txt'), title='DeepGreen-Poisson-Plate2D')
logger.set_names(['Learning Rate', 'TrainLoss', 'TestLoss', 'TrainL(reglr)', 'TrainL(bndry)', 'TrainL(snglr)', 'TrainL(symtry)'])
     
# train and test 
trainloss_curve, testloss_curve, eigenloss_curve, eigenloss_curve_notrain = [], [], [[]], [[]]

# load trained record for resuming training
if args.resume:
    try:
        trainloss_best = 1e10
        model.Xavier_initi() 
        checkpoint_past = torch.load(os.path.join(args.checkpoint, 'model_best.pth.tar'))
        trainloss_best = checkpoint_past['trainloss_best']
        model.load_state_dict(checkpoint_past['state_dict'])
        optimizer.load_state_dict(checkpoint_past['optimizer'])
        schedular.load_state_dict(checkpoint_past['lr_schedular'])
        print('Resume training from the best results of the last training')
        
    except FileNotFoundError:
        print('No training record found, start a new training',"\n")
        
    except KeyError:
        print('Wrongly stored training record!', "\n")
            
    except Exception:
        print('Uncaught error!',"\n")
        
else:
    trainloss_best = 1e10
    model.Xavier_initi()
    

since = time.time()
for epoch in range(num_epochs):
      
    print('Epoch {}/{}'.format(epoch, num_epochs-1), 'with LR = {:.1e}'.format(optimizer.param_groups[0]['lr']))  

    # execute training and testing
    trainloss_epoch, trainloss_reglr_epoch, trainloss_bndry_epoch, trainloss_snglr_epoch, trainloss_symtry_epoch  = train_epoch(epoch, model, optimizer, device)
    testloss_epoch = test_epoch(epoch, model, optimizer, device)
    
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
                            'testloss_epoch': testloss_epoch,
                            'trainloss_best': trainloss_best,
                            'optimizer': optimizer.state_dict(),
                            'lr_schedular': schedular.state_dict(),
                        }, is_best, checkpoint=args.checkpoint)

  
    # save training process to log file
    logger.append([optimizer.param_groups[0]['lr'], trainloss_epoch, testloss_epoch, trainloss_reglr_epoch, trainloss_bndry_epoch, trainloss_snglr_epoch, trainloss_symtry_epoch])
    
    # adjust learning rate according to predefined schedule
    schedular.step()

    # print results
    trainloss_curve.append(trainloss_epoch)
    testloss_curve.append(testloss_epoch)

    print('==> Full-Batch Training Loss = {:.4e}'.format(trainloss_epoch), '  Testing Loss = {:.4e}'.format(testloss_epoch), "\n")
    
logger.close()
time_elapsed = time.time() - since

print('Done in {}'.format(str(datetime.timedelta(seconds=time_elapsed))), '!')
print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################

##############################################################################################
# plot learning curves
fig = plt.figure()
plt.plot(torch.log10(torch.tensor(trainloss_curve)), c = 'red', label = 'training loss' )
plt.title('Learning Curve during Training')
plt.legend(loc = 'upper right')
# plt.show()
fig.savefig(os.path.join(args.image,'TrainCurve.png'))

fig = plt.figure()
plt.plot(torch.log10(torch.tensor(testloss_curve)), c = 'red', label = 'testing loss' )
plt.title('Learning Curve during Testing')
plt.legend(loc = 'upper right')
# plt.show()
fig.savefig(os.path.join(args.image,'TestCurve.png'))
##############################################################################################

##############################################################################################
print('*', '-' * 45, '*')
print('===> loading trained model for inference ...')

# load trained model
checkpoint = torch.load(os.path.join(args.checkpoint, 'model_best.pth.tar'))
model.load_state_dict(checkpoint['state_dict'])

# compute NN predicution for a fixed y
with torch.no_grad():    
    SmpPts_Plot = plotdata.SmpPts_Plot
    SmpPts_Plot = SmpPts_Plot.to(device)
    NN_Plot = model(torch.cat((SmpPts_Plot,Exact_Solution.augmented_func(SmpPts_Plot)),dim=1)) 
    NN_Plot = NN_Plot.to('cpu')

# plot results
x = plotdata.SmpPts_Plot[:,0].reshape(num_test_pts,num_test_pts)
y = plotdata.SmpPts_Plot[:,1].reshape(num_test_pts,num_test_pts)

with torch.no_grad(): 

    # plot exact solution using testing sample points
    plt.figure()
    G_Exact = plotdata.G_Exact_SmpPts.reshape(num_test_pts,num_test_pts)
    G_Exact[torch.isnan(G_Exact)] = 1e3
    G_Exact[G_Exact>1e3] = 1e3
    plt.contourf(x, y, G_Exact, 100, cmap = 'jet')
    plt.title('Exact Green Function on Test Dataset')
    bounds = torch.linspace(G_Exact.min(), G_Exact.max(), 5)
    plt.colorbar(ticks=bounds)
    # plt.show()  
    plt.savefig(os.path.join(args.image,'Green_Exact_TestData.png'))
    plt.close(fig)

    plt.figure()
    G_NN = NN_Plot.reshape(num_test_pts, num_test_pts)  
    G_NN[G_NN > 1e3] = 1e3
    plt.contourf(x, y, G_NN, 40, cmap = 'jet')
    plt.title('Deep Green Function on Test Dataset')
    bounds = torch.linspace(G_NN.min(), G_NN.max(), 5)
    plt.colorbar(ticks=bounds)
    # plt.show()  
    plt.savefig(os.path.join(args.image,'Green_NN_TestData.png'))
    plt.close(fig)

    plt.figure()
    PtErr = G_Exact - G_NN
    plt.contourf(x, y, PtErr, 40, cmap = 'jet')
    plt.title('Pointwise Error on Test Dataset')
    plt.colorbar()
    # plt.show() 
    plt.savefig(os.path.join(args.image,'Green_PtErr_TestData.png'))
    plt.close(fig)        
        
print('*', '-' * 45, '*', "\n", "\n")
##############################################################################################



