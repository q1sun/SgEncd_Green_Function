import numpy as np
import torch

class Conv_1D():
    def __init__(self, model, device, number_of_elements = 2**10, number_of_quadrature_points = 8, left=0., right=1. ):
        super().__init__()
        self.model = model
        self.device = device
        self.number_of_elements = number_of_elements
        self.number_of_quadrature_points = number_of_quadrature_points
        self.left = left
        self.right = right
        
    def conv(self, Quad_Func=None, with_grad = False):
        model = self.model
        device = self.device
        number_of_elements = self.number_of_elements
        number_of_quadrature_points = self.number_of_quadrature_points
        V_mesh, V_Qpts, W_Qpts = self.generate_mesh_quadrature_points()
        # SmpPts_QuadPts = np.column_stack((np.tile(V_mesh, (number_of_elements * number_of_quadrature_points, 1)).reshape(-1),
        #                            np.tile(V_Qpts, (number_of_elements + 1, 1)).reshape(-1)))
        SmpPts_QuadPts = np.column_stack((np.repeat(V_mesh, number_of_elements * number_of_quadrature_points, axis=1).reshape(-1), np.tile(V_Qpts, (number_of_elements + 1, 1)).reshape(-1)))
        SmpPts_QuadPts_symtry = SmpPts_QuadPts[:, [1, 0]]

        SmpPts_QuadPts_Z = np.abs(SmpPts_QuadPts[:, 1] - SmpPts_QuadPts[:, 0])
        SmpPts_QuadPts = np.column_stack((SmpPts_QuadPts, SmpPts_QuadPts_Z))
        SmpPts_QuadPts_symtry = np.column_stack((SmpPts_QuadPts_symtry, SmpPts_QuadPts_Z))   
        
        if not with_grad:
            
            with torch.no_grad():
            
                SmpPts_QuadPts = torch.from_numpy(SmpPts_QuadPts)
                SmpPts_QuadPts_symtry = torch.from_numpy(SmpPts_QuadPts_symtry)
                V_Qpts = torch.from_numpy(V_Qpts)
                W_Qpts = torch.from_numpy(W_Qpts)
                
                SmpPts_QuadPts = SmpPts_QuadPts.to(device)
                SmpPts_QuadPts_symtry = SmpPts_QuadPts_symtry.to(device)
                W_Qpts = W_Qpts.to(device)
                
                NN_Qpts = model(SmpPts_QuadPts)
                NN_Qpts_Symtry = model(SmpPts_QuadPts_symtry)
                DG_Qpts = 0.5 * ( NN_Qpts + NN_Qpts_Symtry ) 
                
        else:
            
            SmpPts_QuadPts = torch.from_numpy(SmpPts_QuadPts)
            SmpPts_QuadPts_symtry = torch.from_numpy(SmpPts_QuadPts_symtry)
            V_Qpts = torch.from_numpy(V_Qpts)
            W_Qpts = torch.from_numpy(W_Qpts)
            
            SmpPts_QuadPts = SmpPts_QuadPts.to(device)
            SmpPts_QuadPts_symtry = SmpPts_QuadPts_symtry.to(device)
            W_Qpts = W_Qpts.to(device)
            
            NN_Qpts = model(SmpPts_QuadPts)
            NN_Qpts_Symtry = model(SmpPts_QuadPts_symtry)
            DG_Qpts = 0.5 * ( NN_Qpts + NN_Qpts_Symtry ) 
                
            
        if Quad_Func.ndim == 1:
            Quad_Func = Quad_Func[:,None]
            num_repeat = 1
        else:
            num_repeat = Quad_Func.shape[1]
        Quad_Func = Quad_Func.repeat( number_of_elements + 1, 1 )
            
        tmp = DG_Qpts * Quad_Func * W_Qpts.repeat( number_of_elements + 1 , 1)
        
        sol = torch.zeros((number_of_elements+1, num_repeat))
        
        for i in range(num_repeat):
            
            sol[:,i] = torch.sum ( torch.reshape( tmp[:,i:i+1], ( number_of_elements + 1,  number_of_elements * number_of_quadrature_points) ) ,dim=1)
        
        return sol
    

        
    def generate_mesh_quadrature_points(self):
        np.set_printoptions(precision=4, suppress=True)
        left = self.left
        right = self.right
        number_of_elements = self.number_of_elements
        number_of_quadrature_points = self.number_of_quadrature_points
        
        # mesh
        V_mesh = np.linspace(left, right, number_of_elements+1)

        # quadrature points
        quadrature_weight_reference_1D, quadrature_point_reference_1D = self.Quadrature_Rule_Reference_1D(number_of_quadrature_points)

        # quadrature points in each element
        V_Qpts = np.zeros(number_of_elements * number_of_quadrature_points)  # coordinates
        W_Qpts = np.zeros(number_of_elements * number_of_quadrature_points)  # weights

        for n in range(number_of_elements):
            # generate quadrature weights and points on local element
            lower_bound = V_mesh[n]
            upper_bound = V_mesh[n+1]
            quadrature_weight_local_1D, quadrature_point_local_1D = self.Quadrature_Rule_Local_1D(quadrature_weight_reference_1D,
                                                                                            quadrature_point_reference_1D,
                                                                                            lower_bound, upper_bound)

            V_Qpts[number_of_quadrature_points * n: number_of_quadrature_points * (n+1)] = quadrature_point_local_1D
            W_Qpts[number_of_quadrature_points * n: number_of_quadrature_points * (n+1)] = quadrature_weight_local_1D
            
        V_mesh = V_mesh[:,None]
        V_Qpts = V_Qpts[:,None]
        W_Qpts = W_Qpts[:,None]
        
        return V_mesh, V_Qpts, W_Qpts
    
        
    def Quadrature_Rule_Local_1D(self, quadrature_weight_reference_1D, quadrature_point_reference_1D, lower_bound, upper_bound):
        np.set_printoptions(precision=4, suppress=True)
        
        # generate (Gauss) quadrature weights and quadrature points on an arbitrary interval [lower_bound,upper_bound]
        quadrature_weight_local_1D = (upper_bound - lower_bound) / 2 * quadrature_weight_reference_1D
        quadrature_point_local_1D = (upper_bound - lower_bound) / 2 * quadrature_point_reference_1D + (upper_bound + lower_bound) / 2
        
        return quadrature_weight_local_1D, quadrature_point_local_1D


    def Quadrature_Rule_Reference_1D(self, number_of_quadrature_points):
        np.set_printoptions(precision=4, suppress=True)
        
        # generate (Gauss) quadrature weights and quadrature points on the reference interval [-1,1]
        if number_of_quadrature_points == 2:
            quadrature_weight_reference_1D = np.array([1, 1])
            quadrature_point_reference_1D = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])
        elif number_of_quadrature_points == 4:
            quadrature_weight_reference_1D = np.array([0.3478548451, 0.3478548451, 0.6521451549, 0.6521451549])
            quadrature_point_reference_1D = np.array([0.8611363116, -0.8611363116, 0.3399810436, -0.3399810436])
        elif number_of_quadrature_points == 8:
            quadrature_weight_reference_1D = np.array([0.1012285363, 0.1012285363, 0.2223810345, 0.2223810345, 0.3137066459,
                                                    0.3137066459, 0.3626837834, 0.3626837834])
            quadrature_point_reference_1D = np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774,
                                                    0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
        else:
            raise ValueError("Invalid number_of_quadrature_points value.")
        
        return quadrature_weight_reference_1D, quadrature_point_reference_1D