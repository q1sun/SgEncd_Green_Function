import torch
import numpy as np

# exact Green's function
def Green_Exact_Disc2D(x, y):
    r = torch.norm( x - y, dim=1) 
    r_x = torch.norm(x, dim=1, keepdim=True)
    x_symtry = x / r_x ** 2
    r_symtry = torch.norm( y - x_symtry, dim=1)
    rho = torch.norm(x, dim=1)
    return - torch.log(r / (r_symtry * rho)) / ( 2 * np.pi)

def augmented_func(x):
    return torch.log( torch.norm((x[:,0:2] - x[:,2:4]), dim=1, keepdim=True))



