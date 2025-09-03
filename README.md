# Singularity-Encoded Green's Function

Green's function provides an inherent connection between theoretical analysis and numerical methods for elliptic partial differential equations, and general absence of its closed-form expression  necessitates surrogate modeling to guide the design of effective solvers. Unfortunately, numerical computation of Green's function remains challenging due to its doubled dimensionality and intrinsic singularity. In this paper, we present a novel singularity-encoded learning approach to resolve these problems in an unsupervised fashion. Our method embeds the Green's function within a one-order higher-dimensional space by encoding its prior estimate as an augmented variable, followed by a neural network parametrization to manage the increased dimensionality. By projecting the trained neural network solution back onto the original domain, our deep surrogate model exploits its spectral bias to accelerate conventional iterative schemes, serving either as a preconditioner or as part of a hybrid solver. The effectiveness of our proposed method is empirically verified through numerical experiments  with two and four dimensional Green's functions, achieving satisfactory resolution of singularities and acceleration of iterative solvers.

## Citation

    @article{sun2025learning,
      title={Learning Singularity-Encoded Green’s Functions with Application to Iterative Methods},
      author={Qi Sun, Shengyan Li, Bowen Zheng, Lili Ju, and Xuejun Xu},
      journal={to be announced},
      year={2025}
    }
