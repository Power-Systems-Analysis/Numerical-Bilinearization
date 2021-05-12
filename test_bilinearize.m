d11 = -1;
d12 = -1;
d21 = 1;

f1 = @(x) d11 * x(1) + d12 * sin(asin(0.5) + x(2));
f2 = @(x) d21 * x(1);
f = {f1; f2};

b11 = 1;
b1 = @(x) b11;
b2 = @(x) 0;
b = {b1; b2};

dt = 1e-4;
t = 0:dt:20.0;
t_size = size(t, 2);
u = 0.5 - 0.5 * cos(t);

tic;
[A, N, B, C] = bilinearize(f, b);
toc

% Nonlinear.
y_nonlin = zeros(t_size, 2);
x = zeros(2, 1);
for i = 1:t_size
    y_nonlin(i, :) = x;
    dx = [
        f{1}(x) + b{1}(x)*u(i);
        f{2}(x) + b{2}(x)*u(i);
    ];
    x = x + dx * dt;
end

% Linear.
A_1 = A(1:2,1:2);
N_1 = N(1:2,1:2);
B_1 = B(1:2,1);
C_1 = C(1:2,1);
a_size = size(A_1, 1);
y_lin = zeros(t_size, 2);
x = zeros(a_size, 1);
for i = 1:t_size
    y_lin(i, :) = x(1:2);
    dx = A_1 * x + N_1 * x * u(i) + B_1 * u(i) + C_1;
    x = x + dx * dt;
end

% Bilinear 2 order.
A_2 = A(1:6,1:6);
N_2 = N(1:6,1:6);
B_2 = B(1:6,1);
C_2 = C(1:6,1);
a_size = size(A_2, 1);
y_bilin_2 = zeros(t_size, 2);
x = zeros(a_size, 1);
for i = 1:t_size
    y_bilin_2(i, :) = x(1:2);
    dx = A_2 * x + N_2 * x * u(i) + B_2 * u(i) + C_2;
    x = x + dx * dt;
end

% Bilinear 3 order.
a_size = size(A, 1);
y_bilin_3 = zeros(t_size, 2);
x = zeros(a_size, 1);
for i = 1:t_size
    y_bilin_3(i, :) = x(1:2);
    dx = A * x + N * x * u(i) + B * u(i) + C;
    x = x + dx * dt;
end


fig = figure;
hold on;
plot(t, y_nonlin(:,1));
plot(t, y_lin(:,1));
plot(t, y_bilin_2(:,1));
plot(t, y_bilin_3(:,1));
legend('nonlin', 'lin', 'bilin 2', 'bilin 3', 'FontSize', 14);
set(findall(fig, 'Type', 'Line'), 'LineWidth', 1.2);
grid on;
title('x_1');

fig = figure;
hold on;
plot(t, y_nonlin(:,2));
plot(t, y_lin(:,2));
plot(t, y_bilin_2(:,2));
plot(t, y_bilin_3(:,2));
legend('nonlin', 'lin', 'bilin 2', 'bilin 3', 'FontSize', 14);
set(findall(fig, 'Type', 'Line'), 'LineWidth', 1.2);
grid on;
title('x_2');
