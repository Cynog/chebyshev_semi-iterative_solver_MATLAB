% matrix-vector multiplication to test the chebyshev solver
function y = matfun_Q(x, N)
    y = [0; x(1:N)] + [2*x(1); 4*x(2:N); 2*x(N+1)] + [x(2:N+1); 0];
end