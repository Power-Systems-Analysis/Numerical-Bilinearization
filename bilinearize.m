% Bilinearization of the system of differential equations up to the third order.
% dx/dt = F(x) + B(x)u => dx/dt = C + Ax + Nxu + Bu
% Fx - column of nonlinear functions of x.
% Bx - control column (now one control signal).
function [A, N, B, C] = bilinearize(Fx, Bx)
    n = size(Fx, 1); % Number of state.
    p = size(Bx, 2); % Number of controls (now 1).
    % Building the matrix A.
    %    | A1  A2   A3  |
    % A =| A20 A21  A22 |
    %    | 0   A30  A31 |
    [A, A0] = matrix_diff(Fx);
    % Building the constant C.
    %     | A0 |
    % C = | 0  |
    %     | 0  |
    C = zeros(n + n^2 + n^3, 1);
    C(1:n, 1) = A0;
    % Building the matrix N.
    %     | B1   B2   B3  |
    % N = | B20  B21  B22 |
    %     | 0    B30  B31 |
    [N, B0] = matrix_diff(Bx);
    % Building the matrix B.
    %     | B0 |
    % B = | 0  |
    %     | 0  |
    B = zeros(n + n^2 + n^3, p);
    B(1:n) = B0;
end

% Construction of matrices of the first, second and third derivatives with
% respect to the column of functions and the formation of a common matrix.
function [A, A0] = matrix_diff(Fx)
    % Increments, for calculating derivatives.
    d = 1e-9;
    % TODO. Consider choosing a step.
    d1 = d / 2;
    d2 = d^(1/2) / 2;
    d3 = d^(1/3) / 2;
    %
    n = size(Fx, 1);
    x = zeros(n, 1);
    A0 = zeros(n, 1);
    A1 = zeros(n, n);
    A2 = zeros(n, n^2);
    A3 = zeros(n, n^3);
    % Loop over equation.
    for k = 1:n
        f = Fx{k};
        % The value of the constant.
        A0(k, 1) = f(x);
        % The derivative is calculated numerically (central derivative).
        % First order partial derivatives.
        for i = 1:n
            x(i) = d1;
            y_p = f(x);
            x(i) = -d1;
            y_m = f(x);
            x(i) = 0;
            A1(k, i) = y_p - y_m;
        end
        % Second order partial derivatives.
        for i = 1:n
            for j = i:n
                % y_pp
                x(i) = x(i) + d2;
                x(j) = x(j) + d2;
                y_pp = f(x);
                x(i) = 0;
                x(j) = 0;
                % y_pm
                x(i) = x(i) + d2;
                x(j) = x(j) - d2;
                y_pm = f(x);
                x(i) = 0;
                x(j) = 0;
                % y_mp
                x(i) = x(i) - d2;
                x(j) = x(j) + d2;
                y_mp = f(x);
                x(i) = 0;
                x(j) = 0;
                % y_mm
                x(i) = x(i) - d2;
                x(j) = x(j) - d2;
                y_mm = f(x);
                x(i) = 0;
                x(j) = 0;
                %
                val = y_pp - y_pm - y_mp + y_mm;
                A2(k, (i - 1) * n + j) = val;
                if i ~= j
                    A2(k, (j - 1) * n + i) = val;
                end
            end
        end
        % Third order partial derivatives.
        for i = 1:n
            for j = i:n
                for t = j:n
                    % y_ppp
                    x(i) = x(i) + d3;
                    x(j) = x(j) + d3;
                    x(t) = x(t) + d3;
                    y_ppp = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_ppm
                    x(i) = x(i) + d3;
                    x(j) = x(j) + d3;
                    x(t) = x(t) - d3;
                    y_ppm = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_pmp
                    x(i) = x(i) + d3;
                    x(j) = x(j) - d3;
                    x(t) = x(t) + d3;
                    y_pmp = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_pmm
                    x(i) = x(i) + d3;
                    x(j) = x(j) - d3;
                    x(t) = x(t) - d3;
                    y_pmm = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_mpp
                    x(i) = x(i) - d3;
                    x(j) = x(j) + d3;
                    x(t) = x(t) + d3;
                    y_mpp = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_mpm
                    x(i) = x(i) - d3;
                    x(j) = x(j) + d3;
                    x(t) = x(t) - d3;
                    y_mpm = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_mmp
                    x(i) = x(i) - d3;
                    x(j) = x(j) - d3;
                    x(t) = x(t) + d3;
                    y_mmp = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    % y_mmm
                    x(i) = x(i) - d3;
                    x(j) = x(j) - d3;
                    x(t) = x(t) - d3;
                    y_mmm = f(x);
                    x(i) = 0;
                    x(j) = 0;
                    x(t) = 0;
                    %
                    val = y_ppp - y_ppm - y_pmp + y_pmm - y_mpp + y_mpm + y_mmp - y_mmm;
                    A3(k, ((i - 1) * n + (j - 1)) * n + t) = val;
                    if i ~= j || i ~= t
                        A3(k, ((i - 1) * n + (t - 1)) * n + j) = val;
                        A3(k, ((j - 1) * n + (i - 1)) * n + t) = val;
                        A3(k, ((j - 1) * n + (t - 1)) * n + i) = val;
                        A3(k, ((t - 1) * n + (j - 1)) * n + i) = val;
                        A3(k, ((t - 1) * n + (i - 1)) * n + j) = val;
                    end
                end
            end
        end
    end
    % Scale the result with respect to the step and add the coefficient of the Taylor series.
    A1 = A1 / d;
    A2 = A2 / (2 * d);
    A3 = A3 / (6 * d);
    I = eye(n, n);
    II = kron(I, I);
    %
    % A20 = A0 * I + I * A0
    A20 = kron(A0, I) + kron(I, A0);
    % A21 = A1 * I + I * A1
    A21 = kron(A1, I) + kron(I, A1);
    % A22 = A2 * I + I * A2
    A22 = kron(A2, I) + kron(I, A2);
    %
    % A30 = A0 * I * I + I * A0 * I + I * I * A0
    % A30 = (A0 * I + I * A0) * I + I * I * A0
    % A30 = A20 * I + I * I * A0
    A30 = kron(A20, I) + kron(II, A0);
    % A31 = A1 * I * I + I * A1 * I + I * I * A1
    % A31 = (A1 * I + I * A1) * I + I * I * A1
    % A31 = A21 * I + I * I * A1
    A31 = kron(A21, I) + kron(II, A1);
    % Building the matrix A.
    %    | A1  A2   A3  |
    % A =| A20 A21  A22 |
    %    | 0   A30  A31 |
    A = zeros(n + n^2 + n^3, n + n^2 + n^3);
    n2 = n + n * n;
    n3 = n2 + n * n * n;
    A(1:n, 1:n) = A1;
    A(1:n, n+1:n2) = A2;
    A(1:n, n2+1:n3) = A3;
    A(n+1:n2, 1:n) = A20;
    A(n+1:n2, n+1:n2) = A21;
    A(n+1:n2, n2+1:n3) = A22;
    A(n2+1:n3, n+1:n2) = A30;
    A(n2+1:n3, n2+1:n3) = A31;
end
