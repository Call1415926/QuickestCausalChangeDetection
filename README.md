# Quickest Causal Change Point Detection by Adaptive Intervention

MATLAB implementation for the paper **"Quickest Causal Change Point Detection by Adaptive Intervention"**

This repository implements sequential change-point detection in linear Gaussian causal graphs. Its main focus is the use of adaptive interventions to amplify an edge-weight change and reduce the expected detection delay (EDD) under a prescribed average run length (ARL) to false alarm.

The code includes the proposed MAX-AI and MULTI-AI procedures, their randomized-intervention and no-intervention variants, two oracle-intervention references, the GGM and JS-WL-CuSum baselines, synthetic experiments, and the two semi-synthetic case studies described in the paper. All detector files called by `main.m`.

## Contents

- [Method overview](#method-overview)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Configuring experiments](#configuring-experiments)
- [Parameters](#parameters)
- [Simulation and evaluation workflow](#simulation-and-evaluation-workflow)
- [Function reference](#function-reference)
- [Outputs](#outputs)
- [Reproducibility notes](#reproducibility-notes)

## Method overview

### Causal model

The observations follow a linear structural equation model (SEM)

\[
X^t = A^t X^t + U^t,
\qquad
U^t \sim \mathcal{N}(\mu^t,\Sigma^t),
\]

where `A` is a weighted adjacency matrix and `Sigma` is diagonal in the main model. The supplied implementation assumes that the nodes are already arranged in a topological order, so `A` is lower triangular and `A(i,j) ~= 0` means that node `j` is a parent of node `i`.

The pre-change model is specified by `(A, mu, Sigma)`, while the post-change model is specified by `(A_oc, mu_oc, Sigma_oc)`. The main experiments consider a single edge-weight change

\[
\widetilde A = A + \Delta E_{k,j},
\]

although two fixed multiple-change examples are also included.

### Intervention and centralization

An intervention `do(X_j = c_j)` removes all incoming edges of node `j`, sets its noise variance to zero, and fixes its value at `c_j`. `CalDoValue.m` constructs the intervention-value vector

\[
C=(c_1,\ldots,c_p)^\mathsf{T}
\]

from the known pre-change model, the minimum design change magnitude, and the intervention gap `delta`. These values are designed so that, under the single-edge-change setting, the origin node of the changed edge has the largest post-change Kullback-Leibler (KL) divergence.

Before monitoring, each observation is centralized and standardized using the pre-change SEM. Under the pre-change model, the resulting variables follow a standard multivariate normal distribution. This transformation is the basis of both proposed detectors.

### Detection and intervention strategies

The core function `getW.m` combines one of two monitoring statistics with one of three intervention policies.

| Code setting | Method label | Description |
| --- | --- | --- |
| `monit_type = 'max'`, `inter_type = 'AI'` | MAX-AI | Node-wise likelihood-ratio CUSUM with adaptive intervention |
| `monit_type = 'max'`, `inter_type = 'RI'` | MAX-RI | Node-wise likelihood-ratio CUSUM with randomized intervention |
| `monit_type = 'max'`, `inter_type = 'NI'` | MAX-NI | Node-wise likelihood-ratio CUSUM without intervention |
| `monit_type = 'multi'`, `inter_type = 'AI'` | MULTI-AI | Joint likelihood-ratio CUSUM with adaptive intervention |
| `monit_type = 'multi'`, `inter_type = 'RI'` | MULTI-RI | Joint likelihood-ratio CUSUM with randomized intervention |
| `monit_type = 'multi'`, `inter_type = 'NI'` | MULTI-NI | Joint likelihood-ratio CUSUM without intervention |

MAX-AI estimates node-wise means and variances and takes the maximum over the node-specific CUSUM statistics. It is computationally lighter and can be used with a shorter estimation window.

MULTI-AI estimates the joint mean vector and covariance matrix under each intervention and uses a multivariate likelihood ratio. It requires a longer window but is intended to retain more information when the single-change concentration property is violated.

At exploration times, AI selects uniformly from no intervention and the `p` node interventions. At exploitation times, it selects the candidate with the largest estimated KL divergence. The parameter `q` controls the number of exploration times in each block of length `w`; all included experiments use `eta = 1`.

The repository also considers:

- **MAX-Oracle (MAX-OI)** and **MULTI-Oracle (MULTI-OI):** always intervene on the known origin node of the changed edge;
- **GGM:** the known pre-change precision-matrix sequential detector of Keshavarz et al. (JMLR, 2020), applied to observational data without intervention;
- **JS-WL-CuSum (JSWL):** a fixed-window CUSUM using a positive-part James-Stein estimator after pre-change standardization.

## Requirements

The code is written in MATLAB. A recent MATLAB release is recommended.

Required products:

- MATLAB;
- Statistics and Machine Learning Toolbox, for `mvnrnd`, `mvnpdf`, and `normpdf`;
- Parallel Computing Toolbox, because `main.m` uses `parpool` and `parfor`.

The experiments can be computationally expensive. The paper settings use `Nrep = 1000`, a monitoring horizon as large as `n = 20000`, separate in-control and out-of-control simulations, multiple ARL targets, and ten methods.

## Repository structure

```text
.
|-- main.m                 Main experiment, calibration, evaluation, and plotting script
|-- getW.m                 MAX/MULTI detectors with AI, RI, or NI policies
|-- getW_maxOra.m          MAX detector with the oracle intervention node
|-- getW_multiOra.m        MULTI detector with the oracle intervention node
|-- getW_GGM.m             GGM baseline
|-- getW_JSWL.m            JS-WL-CuSum baseline
|-- CalDoValue.m           Constructs intervention values from the pre-change SEM
|-- DoOperator.m           Applies do(X_j = c_j) to the SEM parameters
|-- CalDisX.m              Computes the distribution of X under an intervention
|-- CalDisY.m              Computes the centralized post-change distribution
|-- normalKL.m             KL divergence from N(mu,Sigma) to N(0,I)
|-- generateSparseDAG.m    Generates a random lower-triangular linear Gaussian SEM
|-- generate_data.m        Samples observational or interventional SEM data
|-- NECinitial.mat         Ecology case-study graph and parameters
|-- mindInitial.mat        Psychology case-study graph and parameters
|-- tt.txt                 Optional progress log written by detector functions
`-- README.md              This document
```

## Quick start

1. Clone or download the repository and set the repository root as the current MATLAB folder.
2. Open `main.m`.
3. Keep exactly one experiment-configuration block active. The synthetic block is active by default; the ecology, psychology, and multiple-change blocks are commented out.
4. Run:

```matlab
main
```

`main.m` starts a four-worker parallel pool when no pool already exists, generates in-control and out-of-control statistic paths, calibrates a threshold for each target ARL, computes the EDDs, and plots EDD against ARL on a logarithmic horizontal axis.

The default full run evaluates all ten methods listed in `method_names`, including the GGM baseline implemented in `getW_GGM.m`.

For a quick code check before launching the paper-scale experiment, reduce `Nrep` and `n`, and shorten `ARL_set`. Keep `n` larger than the largest target ARL. Small smoke-test settings verify execution but are not intended to reproduce the reported results.

## Configuring experiments

### 1. Synthetic experiment

The first block of `main.m` generates a random lower-triangular SEM. The default configuration is

```matlab
p = 6;
weightRange = [1, 2];
maxInDegree = 2;
muRange = [-1, 1];
SigmaRange = [1/2, 2];
ChangeRange = [0.1, 2];
myChange = 0.1;
delta = 1;
w = 80;
q = 40;
n = 20000;
eta = 1;
Nrep = 1000;
```

The script randomly chooses a strictly lower-triangular entry `(change_row, change_col)` and sets

```matlab
A_oc(change_row, change_col) = ...
    A(change_row, change_col) + myChange;
```

The oracle intervention node is `change_col`, i.e., the origin/parent node of the changed edge. If the selected pre-change entry is zero, this setting represents edge appearance; otherwise it represents a change in an existing edge weight.

`generateSparseDAG.m` constructs the code-level random graph by assigning every non-root node a random number of parents between 1 and `min(maxInDegree,i-1)`, then sampling those parents from earlier nodes. Nonzero weights, noise means, and diagonal noise variances are sampled uniformly from the specified ranges.

### 2. Ecology case study

Comment out the synthetic block and uncomment **Set parameters for Ecology Experiment**. This block loads `NECinitial.mat`, which contains:

- `A`: an `11 x 11` weighted adjacency matrix;
- `mu`: an `11 x 1` noise-mean vector;
- `Sigma`: an `11 x 11` noise covariance matrix;
- `colnames`: the 11 variable names.

Two edge changes are provided in `main.m`:

- Change 1: `A_oc(11,3) = A(11,3) + myChange`, with oracle node 3;
- Change 2: `A_oc(11,9) = A(11,9) + myChange`, with oracle node 9.

Activate only one change definition at a time.

### 3. Psychology case study

Comment out the synthetic block and uncomment **Set parameters for Psychology Experiment**. This block loads `mindInitial.mat`, which contains the analogous `A`, `mu`, `Sigma`, and `colnames` variables for the five-facet mindfulness graph.

Two edge changes are provided:

- Change 1: `A_oc(3,1) = A(3,1) + myChange`, with oracle node 1;
- Change 2: `A_oc(5,4) = A(5,4) + myChange`, with oracle node 4.

Activate only one change definition at a time.

### 4. Multiple-change scenarios

Two fixed multiple-edge-change examples from the paper are included as separate commented blocks. These experiments assess empirical robustness beyond the single-change setting used for the main theory. For both examples, the supplied oracle reference intervenes on node 1, which maximizes the joint KL divergence for the specified model.

## Parameters

| Parameter | Meaning |
| --- | --- |
| `p` | Number of nodes |
| `A`, `mu`, `Sigma` | Pre-change SEM parameters |
| `A_oc`, `mu_oc`, `Sigma_oc` | Post-change SEM parameters (`oc` denotes out of control) |
| `ChangeRange` | Design range `[Delta_min, Delta_max]`; `CalDoValue.m` uses its first element |
| `myChange` | Actual edge-change magnitude used to construct `A_oc` |
| `delta` | Required KL-divergence gap used in the intervention-value design |
| `doValue` | Vector of intervention values returned by `CalDoValue.m` |
| `w` | Sliding-window length for estimating post-change parameters |
| `q` | Exploration budget within each block of length `w` |
| `eta` | Exploration exponent; the implementation selects `ceil(q^eta)` exploration times per block |
| `n` | Maximum monitoring horizon for each Monte Carlo replication |
| `Nrep` | Number of independent Monte Carlo replications |
| `tau` | Change time; `Inf` gives an in-control path and `1` gives the worst-case immediate-change path |
| `b` | Detection threshold; `Inf` requests the full statistic path |
| `ARL_set` | Target in-control average run lengths |
| `epsilon_set` | Allowed absolute calibration error for the corresponding ARL targets |
| `oraNode` | Oracle intervention node; for a changed edge `j -> k`, this is `j` |

The paper recommends shorter windows for MAX-type methods and longer windows for MULTI-type methods because the latter estimate intervention-specific covariance matrices. The current `main.m` applies the same selected `w` and `q` to every method in a run so that the methods are compared under the specified experimental setting.

## Simulation and evaluation workflow

The full workflow in `main.m` is:

1. Define or load the pre-change SEM.
2. Construct the post-change SEM and identify the oracle intervention node.
3. Compute `doValue = CalDoValue(A,mu,Sigma,delta,ChangeRange)`.
4. For each method, simulate `Nrep` in-control paths with `tau = Inf`.
5. Simulate `Nrep` out-of-control paths with `tau = 1`.
6. For each target ARL, use bisection on the stored in-control statistic paths to select a threshold.
7. Apply that threshold to the out-of-control paths and compute the mean alarm time as the EDD.
8. Plot EDD versus target ARL.

The detector functions return complete statistic paths because `main.m` calls them with `b_sim = Inf`. Reusing these paths allows several thresholds to be evaluated without rerunning the Monte Carlo simulation for every ARL target.

If no threshold crossing occurs before time `n`, `main.m` records the alarm time as `n`. Thus, `n` acts as a right-censoring horizon for both estimated ARL and EDD. Choose it sufficiently larger than the largest target ARL and the expected detection delays.

## Function reference

### `getW.m`

Implements the six combinations of MAX/MULTI monitoring and AI/RI/NI intervention. It:

- creates the exploration schedule;
- generates observations sequentially from the appropriate pre- or post-change SEM;
- groups the preceding `w` centralized observations by intervention choice;
- estimates intervention-specific post-change parameters;
- chooses randomized or KL-based adaptive interventions;
- updates node-wise or joint CUSUM statistics;
- returns a `1 x n` overall monitoring-statistic path.

For MAX monitoring, the function internally maintains `p` node-wise paths and returns their pointwise maximum.

### `getW_maxOra.m` and `getW_multiOra.m`

Use the fixed `oraNode` intervention throughout the sequence. They provide oracle performance references rather than implementable detectors, since the changed edge is unknown in practice.

### `getW_JSWL.m`

Generates observational data without intervention, standardizes it with the known pre-change parameters, estimates the post-change mean from the previous `w` observations using global-mean positive-part James-Stein shrinkage, and updates a fixed-window Gaussian plug-in CUSUM statistic. This implementation requires `p >= 3`.

### `getW_GGM.m`

Implements the known pre-change precision-matrix version of the sequential Gaussian graphical model (GGM) detector of Keshavarz et al. (JMLR, 2020). It uses observational data only and does not perform interventions.

The routine first removes the known pre-change observational mean

\[
\mu_X=(I-A)^{-1}\mu
\]

and constructs the known pre-change precision matrix

\[
\Omega_0=(I-A)^\mathsf{T}\Sigma^{-1}(I-A).
\]

For each sliding window, it computes the node-wise conditional quadratic forms, applies the barrier function `f(x) = x - 1 - log(x)`, and standardizes their sum using the digamma/trigamma terms and the normalized precision matrix. The returned value is a `1 x n` monitoring-statistic path. Unlike the MAX, MULTI, and JSWL implementations, this routine does not accumulate its window statistic through a CUSUM recursion.

### `CalDoValue.m`

Implements the intervention-value construction corresponding to Algorithm 1 of the paper. Nodes must be topologically ordered. The routine processes nodes in that order and uses previously constructed ancestor intervention values.

### `DoOperator.m`

Applies a perfect intervention to the model parameters by setting row `j` of `A` to zero, setting `mu(j)` to the intervention value, and setting `Sigma(j,j)` to zero. Passing `j = 0` denotes no intervention.

### `CalDisX.m` and `CalDisY.m`

`CalDisX.m` computes the Gaussian mean and covariance of the observed SEM variables under a selected intervention. `CalDisY.m` computes the mean and covariance after centralization by the pre-change model.

### `normalKL.m`

Computes

\[
D\{\mathcal N(\mu,\Sigma)\,\|\,\mathcal N(0,I)\}
=\frac{1}{2}\left(\mu^\mathsf{T}\mu+\operatorname{tr}(\Sigma)-p-\log\det\Sigma\right).
\]

### `generate_data.m`

Samples `n` rows from the observational or interventional Gaussian SEM. The output has size `n x p`.

## Outputs

After a successful run, the main results remain in the MATLAB workspace:

- `EDD_set_10methods`: mean detection delay for every method and ARL target;
- `threshold_set`: calibrated detection thresholds;
- `achieved_ARL_set`: simulated ARLs achieved by those thresholds;
- `ARL_set`: requested ARL targets;
- `method_names`: row labels for the result matrices.

The script also opens an EDD-versus-ARL figure. It does not automatically save the workspace variables or figure. Save them explicitly if the run must be retained, for example:

```matlab
save('simulation_results.mat', ...
    'EDD_set_10methods', 'threshold_set', ...
    'achieved_ARL_set', 'ARL_set', 'method_names');

exportgraphics(gcf, 'edd_vs_arl.pdf', 'ContentType', 'vector');
```

The detector functions append replication indices to `tt.txt`. This file is only a progress log and is not used in the calculations.

## Reproducibility notes


- To make a run reproducible, set the random-number stream before graph generation and configure parallel random streams appropriately for `parfor`.
- Keep exactly one experiment block active. Otherwise, variables from a later block may overwrite an earlier configuration.
- The adjacency matrix must follow the topological node order expected by `CalDoValue.m` and the ancestor calculations in `getW.m`.
- The main theoretical guarantees apply to a single edge-weight change. The multiple-change blocks are empirical robustness experiments.
- Large experiments can take substantial time and memory. Begin with a small smoke test, then restore the paper-scale settings.

