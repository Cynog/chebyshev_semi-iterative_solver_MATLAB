N = 1000000
tol = 1e-14;
maxit = 50;

% function handle for matrix A
fh_poisson1d = @(x) matfun_poisson1d(x, N);

% Function handle for Richardson iteration
fh_richardson = @(v) v;
fh_jacobi = @(v) 1/4 * v;

% Generate a random normally distributed right hand side
b = randn(N, 1);
normb = norm(b);

% Apply Chebyshev with Richardson iteration for the known eigenvalue bounds
eigmax = 6;
eigmin = 2;
tic;
[x_r, flag_r, relres_r, iter_r, resvec_r] = chebyshev(fh_poisson1d, b, fh_richardson, eigmax, eigmin, tol, maxit);
t_r = toc;
relresvec_r = resvec_r / normb;
flag_r
relres_r
iter_r

% Apply Chebyshev with Jacobi iteration for the known eigenvalue bounds (should converge with same speed as 1.5/0.5 = 3 = 6.0/2.0)
eigmax = 1.5;
eigmin = 0.5;
tic;
[x_j, flag_j, relres_j, iter_j, resvec_j] = chebyshev(fh_poisson1d, b, fh_jacobi, eigmax, eigmin, tol, maxit);
t_j = toc;
relresvec_j = resvec_j / normb;
flag_j
relres_j
iter_j

% Plot the relative residuals
figure();
semilogy(0:length(relresvec_r) - 1, relresvec_r, 0:length(relresvec_j) - 1, relresvec_j);
title("Chebyscheff semi-iterative method", "Interpreter", "latex");
legend(append("Richardson   ", sprintf("t = %6fs", t_r)), append("Jacobi   ", sprintf("t = %6fs", t_j)), "Interpreter", "latex");
xlabel("iterations", "Interpreter", "latex");
ylabel("relative $l_2$-residual", "Interpreter", "latex");