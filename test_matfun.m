% script to verify the correct implementation of matfun matrix-vector-multiplications
N = 5;

% 1d Poisson matrix
dim = N;
A = zeros(dim, dim);

for i = 1:dim
    e_i = zeros(dim, 1);
    e_i(i) = 1;

    for j = 1:dim
        e_j = zeros(dim, 1);
        e_j(j) = 1;

        A(i, j) = e_i' * matfun_poisson1d(e_j, N);
    end

end

A
spy(A)

% 2d Poisson matrix
dim = N*N;
AA = zeros(dim, dim);

for i = 1:dim
    e_i = zeros(dim, 1);
    e_i(i) = 1;

    for j = 1:dim
        e_j = zeros(dim, 1);
        e_j(j) = 1;

        AA(i, j) = e_i' * matfun_poisson2d(e_j, N);
    end

end

AA
spy(AA)
