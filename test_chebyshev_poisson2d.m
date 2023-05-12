% parameters
N = 300;
tol = 1e-12;
maxit = 5000;

% function handle for 2d-Poisson matrix
fh_poisson2d = @(x) matfun_poisson2d(x, N);

% Function handle for Richardson iteration
fh_richardson = @(v) v;
fh_jacobi = @(v)1/4 * v;

% Random normally distributed right hand side
b = randn(N * N, 1);
normb = norm(b);

% Apply Chebyshev with Richardson iteration for the known eigenvalue bounds
eigmax = 4 - 4 * cos(pi * N / (N + 1));
eigmin = 4 - 4 * cos(pi / (N + 1));
tic;
[x_r, flag_r, relres_r, iter_r, resvec_r] = chebyshev(fh_poisson2d, b, fh_richardson, eigmax, eigmin, tol, maxit);
t_r = toc;
relresvec_r = resvec_r / normb;
flag_r
relres_r
iter_r

% Apply Chebyshev with Jacobi iteration for the known eigenvalue bounds (should converge with same speed as 1.5/0.5 = 3 = 6.0/2.0)
eigmax = 1 - cos(pi * N / (N + 1));
eigmin = 1 - cos(pi / (N + 1));
tic;
[x_j, flag_j, relres_j, iter_j, resvec_j] = chebyshev(fh_poisson2d, b, fh_jacobi, eigmax, eigmin, tol, maxit);
t_j = toc;
relresvec_j = resvec_j / normb;
flag_j
relres_j
iter_j

% Plot the relative residuals
figure();
semilogy(0:length(relresvec_r) - 1, relresvec_r, 0:length(relresvec_j) - 1, relresvec_j);
title("Chebyscheff semi-iterative method solving 2d-Poisson equation", "Interpreter", "latex");
legend(append("Richardson   ", sprintf("t = %6fs", t_r)), append("Jacobi   ", sprintf("t = %6fs", t_j)), "Interpreter", "latex");
xlabel("iterations", "Interpreter", "latex");
ylabel("relative $l_2$-residual", "Interpreter", "latex");
