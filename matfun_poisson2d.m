% matrix-vector multiplication for the 2d Poisson equation
function y = matfun_poisson2d(x, N)
    y = zeros(N * N, 1);

    y(1:N) = matfun_poisson1d(x(1:N), N) - x(N + 1:2 * N);

    for i = 2:N - 1
        y((i - 1) * N + 1:i * N) = -x((i - 2) * N + 1:(i - 1) * N) + matfun_poisson1d(x((i - 1) * N + 1:i * N), N) - x(i * N + 1:(i + 1) * N);
    end

    y((N - 1) * N + 1:N * N) =- x((N - 2) * N + 1:(N - 1) * N) + matfun_poisson1d(x((N - 1) * N + 1:N * N), N);
end
