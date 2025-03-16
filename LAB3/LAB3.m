%% LAB3
clear all
close all
clc

p_true = [-0.7647, 0.3012, 0, 32.24, 21.41];
Gp = tf(p_true(3:5), [1 p_true(1:2)], 1);

% noiseless simulation
N = 50;
u = unifrnd(0, 1, N, 1);
y = lsim(Gp, u);

% adding noise
dEta = 1;
dEps = 0.01;
eta = unifrnd(-dEta, dEta, N, 1);
eps = unifrnd(-dEps, dEps, N, 1);
y_tilde = y+eta;
u_tilde = u+eps;

% SparsePOP problem definition
np = 5; a = 3; b = 3;
rows = 2*np+2; cols = np+2*N;
ineq_supp = zeros(rows, cols);

A = zeros(rows, np);
B = zeros(rows, N);
C = zeros(rows, N);

% fill A
i = 3;
for j = 1:np
    A(i, j) = 1;
    A(i+1, j) = 1;
    i = i+2;
end

idx = 1;
for k = a:N
    %fill B
    B = zeros(rows, N);
    i = 1;
    for j = a:-1:1
        B(i, j+idx-1) = 0;
        B(i+1, j+idx-1) = 1;
        i = i+2;
    end

    %fill C
    C = zeros(rows, N);
    i = rows-1;
    for j = 1:b
        C(i, j+idx-1) = 0;
        C(i+1, j+idx-1) = 1;
        i = i-2;
    end

    ineq_supp = [A, B, C];

    %coef vector
    coef = zeros(rows, 1);

    j = 0;
    for i = 1:2:2*a
        coef(i) = y_tilde(k-j);
        coef(i+1) = -1;
        j = j+1;
    end

    j = 0;
    for i = 2*a+1:2:rows
        coef(i) = -u_tilde(k-j);
        coef(i+1) = 1;
        j = j+1;
    end

    ineqPolySys{idx}.noTerms = rows;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.dimVar = cols;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineq_supp;
    ineqPolySys{idx}.coef = coef;

    idx = idx+1;
end

%lower bound
lbd = [-1e10*ones(np, 1); -dEta*ones(N, 1); -dEps*ones(N, 1)];

%upper bound
ubd = [1e10*ones(np, 1); dEta*ones(N, 1); dEps*ones(N, 1)];

pmin = zeros(np, 1);
pmax = zeros(np, 1);

for i = 1:np
    param.relaxOrder = 1;
    param.POPsolver = 'active-set';

    % min
    obj_supp = zeros(1, cols);
    obj_supp(i) = 1;
    obj_coef = 1;
    objPoly.noTerms = 1;
    objPoly.dimVar = cols;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;
    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    % max
    obj_supp = zeros(1, cols);
    obj_supp(i) = 1;
    obj_coef = -1;
    objPoly.noTerms = 1;
    objPoly.dimVar = cols;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;
    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end

pmin, pmax
