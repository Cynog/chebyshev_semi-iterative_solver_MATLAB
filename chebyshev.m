function [y, flag, relres, iter, resvec] = chebyshev(A, b, M, eigmax, eigmin, tol, maxit, x0)
    % assume tolerance will not be reached
    flag = 1;

    % standart values for optional inputs
    if (nargin < 6) || isempty(tol)
        tol = 1e-6;
    end

    if (nargin < 7) || isempty(maxit)
        maxit = 20;
    end

    if (nargin < 8) || isempty(x0)
        x0 = zeros(length(b), 1);
    end

    % check if input already is a function handle
    if isa(A, 'function_handle')
        Afun = A;
    else
        Afun = @(x) A * x;
    end

    if isa(M, 'function_handle')
        Mfun = M;
    else
        Mfun = @(x) M * x;
    end

    % check value of tol
    if tol < eps
        warning("Tolerance tol < eps! Setting tol = eps.");
        tol = eps;
    end

    % l2-norm of right hand side
    normb = norm(b);

    % useful factors for iteration
    om = 4 * ((eigmax + eigmin) / (eigmax - eigmin)) ^ 2;
    omega = 2;
    c = 2 / (eigmax + eigmin);

    % preallocate storage for vector of l2-residuals
    resvec = zeros(maxit + 1, 1);

    % set 0-th iterate to starting vector
    y_old = x0;

    % calculate starting residual
    r = b - Afun(x0);

    % check convergence for starting vector
    normr = norm(r);
    resvec(1) = normr;
    relnormr = normr / normb;

    if relnormr <= tol % convergence
        flag = 0;
        y = x0;
        iter = 0;
        relres = relnormr;
        resvec = resvec(1);
        return;
    end

    % set starting vector as iterate with smallest relative residual
    y_min = x0;
    relres_min = relnormr;
    iter_min = 0;

    % apply stationary method to starting residual
    z = Mfun(r);

    % calculate first iterate
    y = c * z + y_old;

    % check convergece for first iterate
    r = b - Afun(y);
    normr = norm(r);
    resvec(2) = normr;
    relnormr = normr / normb;

    if relnormr <= tol % convergence
        flag = 0;
        relres = relnormr;
        iter = 1;
        resvec = resvec(1:2);
        return;
    end

    % check first iterate for smallest relative residual
    if (relnormr < relres_min)
        y_min = y;
        relres_min = relnormr;
        iter_min = 1;
    end

    % apply iterations until maxit is reached
    for k = 2:maxit
        % apply stationary method to last residual
        z = Mfun(r);

        % update factor omega
        omega = om / (om - omega);

        % calculate new iterate
        tmp = y;
        y = omega * (c * z + y) + (1 - omega) * y_old;
        y_old = tmp;

        % calculate new residual
        r = b - Afun(y);

        % check convergence for new iterate
        resvec(k + 1) = norm(r);
        relnormr = resvec(k + 1) / normb;

        if relnormr <= tol % convergence
            flag = 0;
            relres = relnormr;
            iter = k;
            resvec = resvec(1:k + 1);
            return;
        end

        % check current iterate for smallest relative residual
        if (relnormr < relres_min)
            y_min = y;
            relres_min = relnormr;
            iter_min = k;
        end

    end

    % tol not reached by iteration
    warning("Tolerance not reached by iteration! Providing iterate %d with the smallest relative residual %g.", iter_min, relres_min);
    y = y_min;
    relres = relres_min;
    iter = iter_min;
end
