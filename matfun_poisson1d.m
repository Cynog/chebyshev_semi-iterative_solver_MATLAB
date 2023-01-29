% matrix-vector multiplication for the 1d Poisson equation
function y = matfun_poisson1d(x, N)
    y = [0; -x(1:N-1)] + 4*x + [-x(2:N); 0];
end