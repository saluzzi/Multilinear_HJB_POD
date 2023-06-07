### Multilinear_HJB_POD
A multilinear algorithm to approximate the solution of the Hamilton-Jacobi-Bellman (HJB) equation based on the coupling of Tree Structure Algorithm and multilinear HO-POD-DEIM

## Contents

### Numerical test scripts

* `main_Allen_Cahn.m` Optimal control of the Allen-Cahn (AC) equation.
* `main_transport_HO_POD.m` Optimal control of the Advection-Diffusion (AD) equation via HO-POD.
* `main_transport_POD.m` Optimal control of the Advection-Diffusion (AD) equation via POD.
* `main_burgers.m` Optimal control of the 3D Burger's equation.

### Tree construction scripts

* `full_tree_AC.m` Construction of the tree structure and computation of the HO-POD-DEIM basis for the AC equation.
* `full_tree_red_AC.m` Construction the tree structure driven by the reduced dynamics for the AC equation.
* `prun_tree_statistical.m` Construction of the tree strucutre with the statistical pruning criterion for the AC equation.
* `full_tree_pruned_transport.m`  Construction of the tree with geometrical pruning for the matrix version for the AD equation.
* `tree_pod_pruning.m` Construction of the tree with geometrical pruning for the vector version for the AD equation.
* `full_tree_transport_bilinear.m` Computation of the tree structure driven by the reduced dynamics for the AD equation.
* `tree_bilinear.m` Construction of the tree structure for the vector form for the AD equation.
* `tree_bilinear.m` Construction of the tree structure for the vector form for the AD equation.

### ODE solvers

* `mat_SI.m` Computation of the FOM solution for the AC equation.
* `mat_SI_bil` Computation of the FOM solution for the AD equation.
* `mat_SI_pod_transport.m` Computation of the ROM solution for the AD equation.
* `integrate_and_basis_tensor.m` Integrating the full tensor model on the tree and contstructing the HO-POD-DEIM basis
* `red_tree_tensor.m` Integrating the reduced tensor model on the tree

### Auxiliary

* `refiniment_tree.m` Construction of the refined tree for the statistical pruning
* `value_function.m`  Computation of the value function, optimal trajectory and optimal control.
* `check_mat.m` Auxiliary function for the geometrical pruning for the matrix format.
* `check.m` Auxiliary function for the geometrical pruning for the vector format.

### Tensor Helpers

* `tt_subsubsref.m`
* `sthosvd.m`
* `sylvsolve_tensor.m` Solves the linear tensor equation without vectorization

### Datasets

* `offline_AC.mat`  Offline phase data for the AC equation.
* `final_results_AC.mat` Final results for the AC equation.
