import torch
##############################################################################################
# draw sample points uniformly at random 

def SmpPts_Reglr_Line1D(num_pts_x, num_pts_y):
    ReglrPts_x = torch.rand(num_pts_x,1) 
    ReglrPts_y = torch.rand(num_pts_y,1) 
    ReglrPts_x_Reshape =  ReglrPts_x.repeat(num_pts_y, 1)
    ReglrPts_y_Reshape = torch.reshape( ReglrPts_y.repeat(1, num_pts_x), (num_pts_x*num_pts_y, 1) )

    return torch.cat((ReglrPts_x_Reshape, ReglrPts_y_Reshape),1)


def SmpPts_Bndry_Line1D(num_pts_y): 
    BndryPts_x = torch.tensor([[0], [1]])
    IntrrPts_y = torch.rand(num_pts_y,1)
    BndryPts_x_Reshape = BndryPts_x.repeat(num_pts_y,1)
    IntrrPts_y_Reshape = torch.reshape(IntrrPts_y.repeat(1,2), (2*num_pts_y,1))
    
    return torch.cat((BndryPts_x_Reshape,IntrrPts_y_Reshape),1)


def SmpPts_Snglr_Line1D(num_snglr_pts):
    Snglr_x = torch.rand(num_snglr_pts,1)
    
    return Snglr_x.repeat(1,2)

                  
# draw equidistant sample points from the entire domain [0,1] * [0,1] 
def SmpPts_Test_Square2D(num_test_pts):
    """ num_test_pts = total number of sampling points for model testing for each dimension """

    xs = torch.linspace(0, 1, steps=num_test_pts)
    ys = torch.linspace(0, 1, steps=num_test_pts)
    x, y = torch.meshgrid(xs, ys)

    return torch.squeeze(torch.stack([x.reshape(1,num_test_pts*num_test_pts), y.reshape(1,num_test_pts*num_test_pts)], dim=-1))  


##############################################################################################