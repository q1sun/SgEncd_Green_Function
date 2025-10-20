import torch
import numpy as np
##############################################################################################
# draw sample points uniformly at random 
def SmpPts_Reglr_Disc2D(num_pts_x, num_pts_y):
    
    r_x = torch.sqrt(torch.rand(num_pts_x, 1, 1))
    r_y = torch.sqrt(torch.rand(num_pts_y, 1, 1))
    theta_x = 2 * np.pi * torch.rand(num_pts_x, 1, 1)
    theta_y = 2 * np.pi * torch.rand(num_pts_y, 1, 1)
    ReglrPts_x = torch.cat( ( r_x * torch.cos(theta_x), r_x * torch.sin(theta_x) ), dim=2)
    ReglrPts_y = torch.cat( ( r_y * torch.cos(theta_y), r_y * torch.sin(theta_y) ), dim=2) 
    ReglrPts_x_Reshape = torch.reshape( ReglrPts_x.repeat(1, num_pts_y, 1), (num_pts_x * num_pts_y, 2) )
    ReglrPts_y_Reshape = torch.reshape( ReglrPts_y.repeat(num_pts_x, 1, 1), (num_pts_x * num_pts_y, 2) )

    return torch.cat((ReglrPts_x_Reshape,ReglrPts_y_Reshape), dim=1)

def SmpPts_Bndry_Disc2D(num_pts_x, num_pts_y): 
  
    r_y = torch.sqrt(torch.rand(num_pts_y, 1, 1))
    theta_x = 2 * np.pi * torch.rand(num_pts_x, 1, 1)
    theta_y = 2 * np.pi * torch.rand(num_pts_y, 1, 1)
    BndryPts_x = torch.cat( ( torch.cos(theta_x), torch.sin(theta_x) ), dim=2)
    BndryPts_y = torch.cat( ( r_y * torch.cos(theta_y), r_y * torch.sin(theta_y) ), dim=2)
    BndryPts_x_Reshape = torch.reshape( BndryPts_x.repeat(1, num_pts_y, 1), (num_pts_x * num_pts_y, 2) )
    BndryPts_y_Reshape = torch.reshape( BndryPts_y.repeat(num_pts_x, 1, 1), (num_pts_x * num_pts_y, 2) )
    
    return torch.cat( (BndryPts_x_Reshape, BndryPts_y_Reshape), dim=1)

def SmpPts_Snglr_Disc2D(num_pts_x, num_pts_y, num_circles):

    r_y = torch.sqrt(torch.rand(num_pts_y, 1, 1))
    theta_y = 2 * np.pi * torch.rand(num_pts_y, 1, 1)
    SnglrPts_y = torch.cat( ( r_y * torch.cos(theta_y), r_y * torch.sin(theta_y) ), dim=2)
    SnglrPts_y = torch.reshape( SnglrPts_y.repeat(1, num_circles, 1), (num_pts_y * num_circles, 1, 2) )
    # generate x-pts on concentric circles
    max_eps = 0.1
    eps = max_eps * 0.8 ** (torch.linspace(1, num_circles, num_circles))
    r_x = torch.unsqueeze(torch.reshape( torch.ones(num_pts_y,1,1) * eps,(-1,1) ), 1 )
    r_x = r_x.repeat(1, num_pts_x, 1)
    theta_x = 2 * np.pi * torch.rand(num_pts_y * num_circles, num_pts_x, 1)
    
    SnglrPts_y_Reshape = SnglrPts_y.repeat(1, num_pts_x, 1)
    SnglrPts_x_Reshape = SnglrPts_y_Reshape + r_x * torch.cat( (torch.cos(theta_x), torch.sin(theta_x)), dim=2 )
    
    return torch.cat( (SnglrPts_x_Reshape, SnglrPts_y_Reshape), dim=2)


def SmpPts_Test_Disc2D(num_test_pts):
    # num_test_pts = total number of sampling points for model testing for each dimension 

    r_x = torch.sqrt(torch.rand(num_test_pts * num_test_pts, 1))
    r_y = torch.sqrt(torch.rand(num_test_pts * num_test_pts, 1))
    theta_x = 2 * np.pi * torch.rand(num_test_pts * num_test_pts, 1)
    theta_y = 2 * np.pi * torch.rand(num_test_pts * num_test_pts , 1) 
    
    x = torch.cat( (r_x * torch.cos(theta_x), r_x * torch.sin(theta_x)), dim=1)
    y = torch.cat( (r_y * torch.cos(theta_y), r_y * torch.sin(theta_y)), dim=1)

    return torch.cat( (x,y), dim=1)

# plot green's function for a fixed y = [0,0]
def SmpPts_Plot_Disc2D(num_test_pts, y0_test=torch.zeros(1,2)):

    xs1 = torch.linspace(-1, 1, steps=num_test_pts)
    xs2 = torch.linspace(-1, 1, steps=num_test_pts)
    x1, x2 = torch.meshgrid(xs1, xs2)
    x = torch.squeeze(torch.stack([x1.reshape(1,num_test_pts*num_test_pts), x2.reshape(1,num_test_pts*num_test_pts)], dim=-1))

    return torch.cat((x, y0_test.repeat(x.shape[0],1)), dim=1)



##############################################################################################

