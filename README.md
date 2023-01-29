# chebyshev_semi-iterative_solver_MATLAB

## 1. Overview
This repository contains a MATLAB implementation of the Chebyshev semi-iterative solver for usually sparse, symmetric positive definite matrices with known eigenvalue bounds. The algorithm is oriented on [[Hac]](#1).



## 2. Usage
Let us consider the system of linear equations $A x = b$ for a symmetric positive definite $A \in \mathbb{R}^{n \times n}$. Let us further consider a splitting $A = M - N$ where $M$ is symmetric and $M^{-1} N$ has only real eigenvalues that fulfill $\sigma(M^{-1} N) \in [\alpha, \beta]$ with $\alpha\text{, } \beta \in \mathbb{R}$ and $\alpha < \beta < 1$.
Further define $\kappa = \frac{1 - \alpha}{1 - \beta}$ and $c = \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}$. We have
$$\|y^{(k)} - x\|_A \leq \frac{2c^k}{1 + 2c^k} \|y^{(0)} - x\|_A\text{.}$$
Equivalently we can estimate the eigenvalue bounds of $M^{-1} A$ by $eigmin$ and $eigmax$ and set $\kappa = \frac{eigmax}{eigmin}$. This follows from $M^{-1} N = M^{-1} (M - A) = I - M^{-1} A$. These values $eigmin$ and $eigmax$ are required parameters for the chebyshev algorithm.



## 3. Examples

### 3.1 Poisson Equation
*in progress*

### 3.2 Chebyshev semi-iterative method as preconditioner
*in progress*



## 4. References
<a id="1">[Hac]</a> W. Hackbusch: Iterative Solution of Large Sparse Systems of Equations, Second Edition, Springer, 2016.