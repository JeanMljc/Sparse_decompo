## ⚙️ Configuration

Before running the project, you need to configure the MATLAB environment by editing the `config.m` script.

This script defines and adds the required folders and solver paths to the MATLAB search path.

To run this code you need **CPLEX 12.10** and **Matlab r2019b**

---

| Variable | Description |
|-----------|--------------|
| `folderPath` | Path to the main project folder to add all subfolders |
| `CPLEXPath` | Path to IBM **CPLEX 12.10** |
| `ElvanderPath` | Path to the code from [Elvander20] |

---

[Elvander20] : https://www.sciencedirect.com/science/article/pii/S0165168420300207 

code : https://github.com/filipelv/felvander.github.io/raw/refs/heads/main/multimarginal_omt.zip


## Scripts

Each script follows the same overall workflow:

1. **Generate ground-truth multichannel data**  
   Build ground-truth data Z*

2. **Add noise to the ground-truth data**  
   Create the noisy observations Z from Z*.

3. **Construct the dictionary**  
   Build the dictionary matrix A (convolutive)

4. **Generate ground-truth signal and decomposition**  
   Define ground-truth decomposition in the dictionary : X_GT

5. **Set solver and model hyperparameters**  

6. **Run the optimization algorithm**  

7. **Plot decompositions $\hat{X}$**  

### Files \& description

|File | description |
|------------------|-------------|
| `config.m` | Configuration file **to complete with your own paths**.|
| `Main_deconv.m` | Solves a multichannel deconvolution problem using proposed approaches (convex and non-convex) : 1ST & 2ND order approaches + lp-norm vs. LS + lp-norm.|
| `Main_deconv_Elvander.m` | Solves a multichannel deconvolution problem using the method proposed in [Elvander20] and the associated code. |
| `Main_deconv_Elvander_tests.m` | Compare several results for multichannel deconvolution problem with [Elvander20] method for different hyperparameters. |
| `R2_interpolation.m` | Solve R2 interpolation problem numerically **with/without overcomplete** convolutive dictionary (oversampled grid) |
| `R2_cost_discrete_line.m` | Compute R2 cost associated to a discrete line $f(\alpha)$. Plot the cost w.r.t to the slop $\alpha$. |
| `Plot_arrow_OT.m` |  to visualise the mass displacement. |
| `Plot_extremals_solutions_interpo_R2.m` | Plot extremal solutions of R2 interpolation problem. |

## 📂 Folder `/Build_data/`

* Build data for the multichannel deconvolution problem : linear and quadratic trajectories, balanced data.  


## 📂 Folder `/Algos/`

### Algorithms for convex criteria

* `Algo_CPLEX_1ST_QP.m` solves convex QP problem : 

$$ Q_1(x,p,w) := \min_{(x,p) \in \mathcal{C}_1} \quad J(x) + \alpha L_1(p) + \alpha' w^T x  $$

* `Algo_CPLEX_2ND_QP.m` solves convex QP problem  : 

$$ Q_2(x,p,w) := \min_{(x,p) \in \mathcal{C}_2} \quad J(x) + \alpha L_2(p) + \alpha' w^T x  $$

* `Algo_CPLEX_WL1.m` solves _weighted_ $\ell_1$ problem with positivity constraints : 

$$ \min_{x \geq 0} \quad J(x) + \lambda w^T x $$

With uniforms weights i.e. $w = 1$, this problem is equivalent to BPDN/LASSO with positivity constraints.

These three algorithms use `cplexqp()` function from **CPLEX 12.10** 

More details at : 

* Annnexe C [These Mallejac] 
* GRETSI25 (https://gretsi.fr/data/colloque/pdf/2025_soussen1567.pdf)


### Algorithms for non-convex criteria 

* `IRL1_1ST_cplex.m` solves non-convex problem : 

$$\min_{(x,p) \in \mathcal{C}_1} \quad J(x) + \alpha L_1(p) + \alpha' ||x||_p$$

* `IRL1_2ND_cplex.m` solves non-convex problem : 

$$\min_{(x,p) \in \mathcal{C}_1} \quad J(x) + \alpha L_2(p) + \alpha' ||x||_p$$

* `IRL1_LP_cplex.m` solves non-convex problem :

$$\min_{x \geq 0} \quad J(x) + \lambda ||x||_p$$

These functions are implementations of IRL1 algorithmic scheme : each iteration amounts to solving $Q_1(x,p,w), Q_2(x,p,w)$ or _weighted_ $\ell_1$ problem under positivity constraints. We use `Algo_CPLEX_1ST_QP.m`, `Algo_CPLEX_2ND_QP.m`, `Algo_CPLEX_WL1.m` to do so.

#### Other approaches

| Function | Description |
|--------------------|-------------|
| `Algo_static_Elvander.m`| Code associated to static method in [Elvander20] |
| `Algo_dynamic_Elvander.m`| Code associated to dynamic method in [Elvander20] |

#### 📂 Folder `Algos/Algo_Janati/`

* `Algo_Janati_BCD.m` solves the problem :
    
    **to be completed** 

    with the algorithm scheme proposed in [Janati19]. The latter proposed to optimize alternatively on inverse problem variables $x_1 ... x_N$ and OT variables $P_1 ... P_{N-1}$

| Function | Description |
|--------------------|-------------|
| `descent_prox.m` | Proximal descent algorithm |
| `OT_sinkhorn_par.m` | Sinkhorn algorithm |
| `sum_MV.m` | Sum function |

---
# Sparse_decompo
