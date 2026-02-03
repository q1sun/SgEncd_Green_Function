# Singularity-Encoded Green's Function

Green's function provides an inherent connection between theoretical analysis and numerical methods for elliptic partial differential equations, and general absence of its closed-form expression  necessitates surrogate modeling to guide the design of effective solvers. Unfortunately, numerical computation of Green's function remains challenging due to its doubled dimensionality and intrinsic singularity. 

In this paper, we present a novel singularity-encoded learning approach to resolve these problems in an unsupervised fashion. Our method embeds the Green's function within a one-order higher-dimensional space by encoding its prior estimate as an augmented variable, followed by a neural network parametrization to manage the increased dimensionality. By projecting the trained neural network solution back onto the original domain, our deep surrogate model exploits its spectral bias to accelerate conventional iterative schemes, serving either as a preconditioner or as part of a hybrid solver. The effectiveness of our proposed method is empirically verified through numerical experiments with two and four dimensional Green's functions, achieving satisfactory resolution of singularities and acceleration of iterative solvers.

---

## Overview

This repository contains the implementation of **singularity-encoded Green’s functions** for partial differential equation (PDE) solvers.  
It includes both the **training framework** for learning the Green’s function and the **downstream applications** that employ the trained Green’s function to construct neural preconditioners, hybrid iterative solvers, and multigrid accelerations.

This code was developed and maintained by Shengyan Li. For any questions or requests, please contact **2410285@tongji.edu.cn**.

It is organized into two main components:

1. **`SgEncd_Green_<Experiment>`**
   
   This module contains scripts for training and evaluating the singularity-encoded Green’s function.
   - Training script: 
     `Train_SgEncd_Green_<Experiment>/train_SgEncd_Green_<Experiment>.py`
   - Testing and visualization (matlab):
     `run_DownStreamTask_SgEncd_Green_<Experiment>.m`
   

2. **`Acltd_IteM_<Experiment>`**
   - This module contains downstream applications using the trained Green’s function to accelerate iterative solvers:
     | Script | Purpose |
     |---------|----------|
     | `run_Neural_Preconditioner_<Experiment>.m` | Neural preconditioning of linear systems |
     | `run_Hybrid_IteM_<Experiment>.m` | Hybrid iterative method combining neural and classical solvers |
     | `run_Hybrid_Multigrid_<Experiment>.m` | Hybrid scheme integrated into multigrid framework |

   
For reproducibility, we summarize the key training hyper-parameters used for learning the singularity-encoded Green’s functions in all experiments.
### Training hyper-parameters for learning Green’s functions

|Two-dimensional<br>Green's Function| Training<br>Epochs | Network<br>(Depth, Width) | Penalty Coefficients<br>($\beta_{\text{Snglr}}, \beta_{\text{Bndry}}, \beta_{\text{Symtr}}$) | Mini-batch<br>Size | Learning Rate<br>(Init., Decay) | AdamW<br>($\beta_1, \beta_2$) |
|:----------------------:|:-----------:|:------------------:|:-------------------------:|:-------------------------------------------------------------------------------------------:|:------------------:|:-------------------------------:|
| Poisson Eq. (15) | 20k | (4, 40) | (400, 400, 400) | 1,600 | ($10^{-3}$, 0.1) | (0.9,&nbsp;0.999) |
| Helmholtz Eq. (33) | 50k | (4, 40) | (400, 400, 400) | 12,500 | ($10^{-3}$, 0.1) | (0.9,&nbsp;0.999) |
| Convection–Diffusion Eq. (40) | 35k | (4, 40) | (400, 400, 400) | 7,500 | ($10^{-3}$, 0.1) | (0.9,&nbsp;0.999) |
|**Four-dimensional<br>Green's Function**| | | | | | |
| Poisson Eq. (36) | 30k | (6, 40) | (400, 400, 400) | 12,800 | ($10^{-3}$, 0.1) | (0.9,&nbsp;0.999) |
| Elliptic Eq. (44) | 30k | (6, 40) | (400, 400, 400) | 12,800 | ($10^{-3}$, 0.1) | (0.9,&nbsp;0.999) |




## Citation

    @article{sun2025learning,
      title={Learning Singularity-Encoded Green’s Functions with Application to Iterative Methods},
      author={Qi Sun, Shengyan Li, Bowen Zheng, Lili Ju, and Xuejun Xu},
      journal={arXiv preprint arXiv:2509.11580},
      year={2025}
    }
