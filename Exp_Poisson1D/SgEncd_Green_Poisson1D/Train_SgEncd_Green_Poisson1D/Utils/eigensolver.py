import numpy as np
import torch
from scipy.linalg import eigh
from scipy import io
import os

class eigensolver_Green_Poisson1D():
    def __init__(self, eigen_N, model, device, basis_type=101, h = 1./2**8, left=0., right=1. ):
        super().__init__()
        self.eigen_N = eigen_N
        self.basis_type = basis_type
        self.h = h
        self.left = left
        self.right = right
        self.model = model
        self.device = device

        
        
        
    def eigensolver(self):
        np.set_printoptions(precision=6, suppress=True)
        
        P_mesh, T_mesh, V_basis, T_basis = self.Generate_Information_Matrices_Linear_Quadratic_1D()
        SmpPts_FEM_grids = np.column_stack((np.tile(V_basis, (V_basis.shape[1], 1)).reshape(-1), np.repeat(V_basis.T, V_basis.shape[1], axis=1).reshape(-1)))
        SmpPts_FEM_grids_Symtry = SmpPts_FEM_grids
        SmpPts_FEM_grids = SmpPts_FEM_grids[:,[1, 0]]
        io.savemat(os.path.join('..','Task2_Represent_Solution','data.mat'),{'SmpPts_FEM_grids':SmpPts_FEM_grids})
        # 1. FEM Setting
        #####################################################################
        # 1-1. generate information and mesh matrices
        #####################################################################
        number_of_elements =T_mesh.shape[1]
        number_of_nodes = V_basis.shape[1]
        
        if self.basis_type == 101:  # linear elements
            number_of_basis = 2
        elif self.basis_type == 102:  # quadratic elements
            number_of_basis = 3
        else:
            raise ValueError("Invalid basis_type value.")
        
        #####################################################################
        # 1-2. quadrature's points and weights on reference interval
        #####################################################################
        number_of_quadrature_points = 8
        quadrature_weight_reference_1D, quadrature_point_reference_1D = self.Quadrature_Rule_Reference_1D(number_of_quadrature_points)
        
        #####################################################################

        # (a) solve generalized eigenvalue problem of KLE
        #####################################################################
        # (a-1) assemble stiffness matrices
        #####################################################################
        lBaG = np.zeros((number_of_basis, number_of_quadrature_points, number_of_elements))
        H = np.zeros((number_of_basis, number_of_basis, number_of_elements))
        K = np.zeros((number_of_nodes, number_of_nodes))
        
        # stiffness matrix
        for n in range(number_of_elements):
            # generate quadrature weights and points on local element
            vertices = P_mesh[:, T_mesh[:, n]-1].flatten()
            lower_bound = min(vertices[0], vertices[1])
            upper_bound = max(vertices[0], vertices[1])
            quadrature_weight_local_1D, quadrature_point_local_1D = self.Quadrature_Rule_Local_1D(quadrature_weight_reference_1D,
                                                                                            quadrature_point_reference_1D,
                                                                                            lower_bound, upper_bound)
            
            for i in range(number_of_basis):
                for j in range(number_of_quadrature_points):
                    # lBaG[i, j, n] denotes the value corresponding to ith local base of nth element at local jth Gauss point
                    lBaG[i, j, n] = self.Local_Basis_1D(quadrature_point_local_1D[j], vertices, self.basis_type, i+1, 0)
            
            for k in range(number_of_basis):
                H[k, :, n] = quadrature_weight_local_1D @ (np.multiply(np.repeat(lBaG[k:k+1, :, n],number_of_basis,axis=0), lBaG[:, :, n])).T
            K[np.ix_(T_basis[:, n]-1, T_basis[:, n]-1)] = H[:, :, n] + K[np.ix_(T_basis[:, n]-1, T_basis[:, n]-1)]
            
            
        with torch.no_grad():
            SmpPts_FEM_grids = torch.from_numpy(SmpPts_FEM_grids)
            SmpPts_FEM_grids_Symtry = torch.from_numpy(SmpPts_FEM_grids_Symtry)
            
            SmpPts_FEM_grids_aug = torch.cat((SmpPts_FEM_grids, torch.abs(SmpPts_FEM_grids[:,0:1]-SmpPts_FEM_grids[:,1:2])), dim=1)
            SmpPts_FEM_grids_Symtry_aug = torch.cat((SmpPts_FEM_grids_Symtry, torch.abs(SmpPts_FEM_grids_Symtry[:,0:1]-SmpPts_FEM_grids_Symtry[:,1:2])), dim=1)
            
            SmpPts_FEM_grids_aug = SmpPts_FEM_grids_aug.to(self.device)
            NN_FEM_grids = self.model(SmpPts_FEM_grids_aug) 

            SmpPts_FEM_grids_Symtry_aug = SmpPts_FEM_grids_Symtry_aug.to(self.device)
            NN_FEM_grids_Symtry = self.model(SmpPts_FEM_grids_Symtry_aug) 
            
            Matrix_DG = ( 0.5 * (NN_FEM_grids + NN_FEM_grids_Symtry) ).cpu().numpy()
            
        # io.savemat(os.path.join('..','Task2_Represent_Solution','data.mat'),{'SmpPts_FEM_grids_aug':SmpPts_FEM_grids_aug.cpu().numpy()})
        
        Matrix_DG = np.reshape(Matrix_DG, (number_of_nodes, number_of_nodes), order='F')
        M = K @ Matrix_DG @ K
        io.savemat(os.path.join('..','Task2_Represent_Solution','data.mat'),{'K':K,'Matrix_DG':Matrix_DG})
        
        #####################################################################
        # (a-2) solve generalized eigenvalue problem
        #####################################################################
        eigenvalues, eigenvectors = eigh(M, K, subset_by_index=[number_of_nodes - self.eigen_N, number_of_nodes - 1])
        mu = np.sort(np.real(eigenvalues))[::-1]
        phi = np.real(eigenvectors[:, np.argsort(np.real(eigenvalues))[::-1]])
        
        # normalize eigenfunctions
        normalize_matrices = np.linalg.norm(phi, axis=0)
        normalize_matrices = np.tile(normalize_matrices, (number_of_nodes, 1))
        phi /= normalize_matrices
        
        #####################################################################
        
        return mu.copy(), phi.copy()


    def Local_Basis_1D(self, x, vertices, basis_type, basis_index, derivative_degree):
        np.set_printoptions(precision=4, suppress=True)
        
        if basis_type == 101:  # linear element
            if derivative_degree == 0:
                if basis_index == 1:
                    result = (vertices[1] - x) / (vertices[1] - vertices[0])
                elif basis_index == 2:
                    result = (x - vertices[0]) / (vertices[1] - vertices[0])
            elif derivative_degree == 1:
                if basis_index == 1:
                    result = 1 / (vertices[0] - vertices[1])
                elif basis_index == 2:
                    result = 1 / (vertices[1] - vertices[0])
        
        elif basis_type == 102:  # quadratic element
            denominator = (vertices[0] - vertices[1])**2
            
            if derivative_degree == 0:
                if basis_index == 1:
                    result = 2*x**2 - (vertices[0] + 3*vertices[1])*x + vertices[1]**2 + vertices[0]*vertices[1]
                elif basis_index == 2:
                    result = 2*x**2 - (3*vertices[0] + vertices[1])*x + vertices[0]**2 + vertices[0]*vertices[1]
                elif basis_index == 3:
                    result = -4*x**2 + 4*(vertices[0] + vertices[1])*x - 4*vertices[0]*vertices[1]
            elif derivative_degree == 1:
                if basis_index == 1:
                    result = 4*x - (vertices[0] + 3*vertices[1])
                elif basis_index == 2:
                    result = 4*x - (3*vertices[0] + vertices[1])
                elif basis_index == 3:
                    result = -8*x + 4*(vertices[0] + vertices[1])
            elif derivative_degree == 2:
                if basis_index == 1:
                    result = 4
                elif basis_index == 2:
                    result = 4
                elif basis_index == 3:
                    result = -8
            
            result = result / denominator
        
        return result




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


    def Generate_Information_Matrices_Linear_Quadratic_1D(self):
        np.set_printoptions(precision=4, suppress=True)
        left = self.left
        right = self.right
        h = self.h
        basis_type = self.basis_type
        
        # 1. generate information matrices P and T for a partition of [left, right]
        ####################################################################
        # 1-1. linear elements (also the grid points)
        ####################################################################
        N_linear = int((right - left) / h)
        P_linear = np.linspace(left, right, N_linear+1)
        T_linear = np.column_stack((np.arange(1, N_linear+1), np.arange(2, N_linear+2)))

        # boundary nodes (global index among all finite element nodes)
        BN_linear = [1, len(P_linear)]
        ####################################################################

        # 1-2. quadratic elements
        ####################################################################
        N_quadratic = int((right - left) / h)
        dh = h / 2
        dN = N_quadratic * 2
        P_quadratic = np.linspace(left, right, dN+1)
        T_quadratic = np.column_stack((np.arange(1, 2*dN, 2), np.arange(3, 2*dN+2, 2), np.arange(2, 2*dN+1, 2)))

        # boundary nodes (global index among all finite element nodes)
        BN_quadratic = [1, len(P_quadratic)]
        ####################################################################

        ####################################################################
        # mesh matrices are the same as the linear elements
        P_mesh = P_linear[:,None]
        T_mesh = T_linear

        # finite elements
        if basis_type == 101:  # linear elements
            V_basis = P_linear[:,None]
            T_basis = T_linear
        elif basis_type == 102:  # quadratic elements
            V_basis = P_quadratic[:,None]
            T_basis = T_quadratic
        else:
            raise ValueError("Invalid basis_type value.")
        ####################################################################

        return P_mesh.transpose(), T_mesh.transpose(), V_basis.transpose(), T_basis.transpose()


