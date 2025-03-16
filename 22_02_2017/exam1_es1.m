clear all 
close all
clc

load('data_exam_2A.mat');

u = Input;
y = Output;

dEta = 0.05;
dEps = 0.05;

N = 30;

%SparsePOP problem definition

% 1
np = 3; dimVar = np+2*N; noTerms1 = 8;
ineqSupp1 = zeros(noTerms1, dimVar);
A = zeros(noTerms1, np); 
B = zeros(noTerms1, N);
C = zeros(noTerms1, N);

%fill time invariant part
A(2:4, 1:3) = eye(3);
A(6:8, 1:3) = eye(3);

idx = 1;
for k = 2:N
    B = zeros(noTerms1, N);
    B(5, k) = 1;
    B(6, k-1) = 1;

    C = zeros(noTerms1, N);
    C(7, k) = 1;
    C(8, k-1) = 1;

    ineqSupp1 = [A, B, C];

    coef1 = zeros(noTerms1, 1);
    coef1 = [y(k:-1:k-1); -u(k:-1:k-1); -1; -1; 1; 1];
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.noTerms = noTerms1;
    ineqPolySys{idx}.supports = ineqSupp1;
    ineqPolySys{idx}.coef = coef1;
    idx = idx+1;
end

%2

noTerms2 = 2;
ineqSupp2 = zeros(noTerms2, dimVar);
ineqSupp2(1:2, 2:3) = 2*eye(2);

coef2 = [1; -1];

ineqPolySys{idx}.typeCone = 1;
ineqPolySys{idx}.degree = 2;
ineqPolySys{idx}.dimVar = dimVar;
ineqPolySys{idx}.noTerms = noTerms2;
ineqPolySys{idx}.supports = ineqSupp2;
ineqPolySys{idx}.coef = coef2;
idx = idx+1;

%3
noTerms3 = 4;
ineqSupp3 = zeros(noTerms3, dimVar);
ineqSupp3(2:4, 1:3) = eye(3);

coef3 = [-3; -3; 1; 1];

ineqPolySys{idx}.typeCone = 1;
ineqPolySys{idx}.degree = 1;
ineqPolySys{idx}.dimVar = dimVar;
ineqPolySys{idx}.noTerms = noTerms3;
ineqPolySys{idx}.supports = ineqSupp3;
ineqPolySys{idx}.coef = coef3;

lbd = [-1; -1e10*ones(np-1, 1); -dEta*ones(N, 1); -dEps*ones(N, 1)];
ubd = [1; 1e10*ones(np-1, 1); dEta*ones(N, 1); dEps*ones(N, 1)];

pmin = zeros(np, 1);
pmax = zeros(np, 1);
for i = 1:np
    param.POPsolver = 'active-set';
    param.relaxOrder = 1;

    %min
    objPoly.noTerms = 1;
    objPoly.dimVar = dimVar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = 1; 
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    %max
    objPoly.noTerms = 1;
    objPoly.dimVar = dimVar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = -1; 
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end

[pmin, pmax]



