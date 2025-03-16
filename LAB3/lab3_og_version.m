%% LAB3
clear all
close all
clc

p_true = [-0.7647, 0.3012, 0, 32.24, 21.41];
Gp = tf(p_true(3:5), [1 p_true(1:2)], 1);
%% noiseless simulation
N = 50;
u = unifrnd(0, 1, N, 1);
y = lsim(Gp, u);

%% adding noise
dEta = 1;
dEps = 0.01;
eta = unifrnd(-dEta, dEta, N, 1);
eps = unifrnd(-dEps, dEps, N, 1);
y_tilde = y + eta;
u_tilde = u + eps;

%% SparsePOP problem definition
params = 5;
a = 3; % eta samples
b = 3; % epsilon samples

% define inequality support matrix
rows = 2*params+2; cols = params+2*N; 
ineq_supp = zeros(rows, cols);

% it can be divided in 3 matrices (theta | eta | epsilon)
A = zeros(rows, params);
B = zeros(rows, N);
C = zeros(rows, N);

% fill A (same for all k)
i = 3;
for j = 1:params
    A(i:i+1, j) = [1; 1];
    i = i+2;
end

idx = 1;
for k = a:N
    % fill B
    B = zeros(rows, N);
    i = 1;
    
    for j = a:-1:1
        B(i:i+1, j+idx-1) = [0; 1];
        i = i+2;
    end
    % fill C
    C = zeros(rows, N);
    i = rows-1;
    for j = 1:1:b
        C(i:i+1, j+idx-1) = [0; 1];
        i = i-2;
    end
    
    ineq_supp = [A, B, C];
    
    % define coeff vector
    coef = zeros(rows, 1);
    
    % fill coeff vector
    j = 0;
    for i = 1:2:2*a
        coef(i) = y_tilde( k-j );
        coef(i+1) = -1;
        j = j+1;
    end
    
    j = 0;
    for i = (2*a+1):2:rows
        coef(i) = -u_tilde( k-j );
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
lbd = [-1e10*ones(params, 1); -dEta*ones(N, 1); -dEps*ones(N, 1)];
%upper bound
upd = [1e10*ones(params, 1); dEta*ones(N, 1); dEps*ones(N, 1)];

p_min_rel = zeros(params, 1);
p_min_ref = zeros(params, 1);
p_max_rel = zeros(params, 1);
p_max_ref = zeros(params, 1);
for i = 1:5
    
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
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, upd, param);
    p_min_rel(i) = POP.xVect(i);
    p_min_ref(i) = POP.xVectL(i);

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
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, upd, param);
    p_max_rel(i) = POP.xVect(i);
    p_max_ref(i) = POP.xVectL(i);
end

[p_min_rel, p_max_rel]
[p_min_ref, p_max_ref]













