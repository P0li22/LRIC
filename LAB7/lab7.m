%% LAB7
clear all
close all
clc

% problem setup
dEta = 0.02;
load('data_exam_2A_hammer.mat');
u_tilde = Input;
y_tilde = Output;
N = 20;

% SparsePOP problem definition
np = 5; nvar = np+2*N; 
ineq1_supp = zeros(6, nvar); % L
ineq2_supp = zeros(3, nvar); % N
ineq3_supp = zeros(4, nvar); % dcgain

A1 = zeros(6, 5); B1 = zeros(6, N); C1 = zeros(6, N);

% fill time invariant parts
A1(3:4, 1) = [1; 1];
A1(5:6, 2:3) = [1 0; 0 1];
ineq2_supp(2:3, 4:5) = [1 0; 0 1];
ineq3_supp(2:4, 1:3) = eye(3);

% fill time variant parts of ineq1 and define inequalities related to L
idx = 1;
for k = 2:N
    %fill B1
    B1 = zeros(6, N);
    B1(6, idx) = 1;
    B1(5, idx+1) = 1;

    % fill C1
    C1 = zeros(6, N);
    C1(2, idx+1) = 1;
    C1(4, idx) = 1;

    ineq1_supp = [A1, B1, C1];

    % coef vector
    coef = [y_tilde(k); -1; y_tilde(k-1); -1; -1; -1];

    ineqPolySys{idx}.noTerms = 6;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.dimVar = nvar;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineq1_supp;
    ineqPolySys{idx}.coef = coef;

    idx = idx+1;
end

%fill time variant parts of ineq2 and define inequalities related to N
for k= 1:N
   %fill first line of ineq2_supp
   ineq2_supp(1, :) = zeros(1, nvar);
   ineq2_supp(1, np+k) = 1;

   coef = [1; -u_tilde(k); -u_tilde(k)^3];

   ineqPolySys{idx}.noTerms = 3;
   ineqPolySys{idx}.degree = 1;
   ineqPolySys{idx}.dimVar = nvar;
   ineqPolySys{idx}.typeCone = -1;
   ineqPolySys{idx}.supports = ineq2_supp;
   ineqPolySys{idx}.coef = coef;

   idx = idx+1;
end

coef = [-1; -1; 1; 1];
ineqPolySys{idx}.noTerms = 4;
ineqPolySys{idx}.degree = 1;
ineqPolySys{idx}.dimVar = nvar;
ineqPolySys{idx}.typeCone = -1;
ineqPolySys{idx}.supports = ineq3_supp;
ineqPolySys{idx}.coef = coef;

%lower bound
lbd = [-1e10*ones(np+N, 1); -dEta*ones(N, 1)];

%upper bound
ubd = [1e10*ones(np+N, 1); dEta*ones(N, 1)];

pmin = zeros(np, 1);
pmax = zeros(np, 1);

for i=1:np
    param.relaxOrder = 2;
    param.POPsolver = 'active-set';

    % min
    obj_supp = zeros(1, nvar);
    obj_supp(i) = 1;
    obj_coef = 1;
    objPoly.noTerms = 1;
    objPoly.dimVar = nvar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;
    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    % max
    obj_supp = zeros(1, nvar);
    obj_supp(i) = 1;
    obj_coef = -1;
    objPoly.noTerms = 1;
    objPoly.dimVar = nvar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;
    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end

pmin, pmax
