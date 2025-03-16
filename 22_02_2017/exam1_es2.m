clear all
close all
clc

load('data_exam_2A_hammer_ver2.mat');
u = Input;
y = Output;
dEta = 0.02;
N = 20;

%sparsePOP problem definition

%1
np = 5; dimVar = np+2*N; noTerms1 = 6;
ineqSupp1 = zeros(noTerms1, dimVar);
A = zeros(noTerms1, np);
B = zeros(noTerms1, N);
C = zeros(noTerms1, N);

% fill time invariant part of A
A(2:4, 1:3) = eye(3);
A(6, 1) = 1;

idx = 1;
for k = 2:N
    B = zeros(noTerms1, N);
    B(3, k) = 1;
    B(4, k-1) = 1;

    C = zeros(noTerms1, N);
    C(5, k) = 1;
    C(6, k-1) = 1;

    ineqSupp1 = [A, B, C];
    
    coef1 = zeros(noTerms1, 1);
    coef1 = [y(k); y(k-1); -1; -1; -1; -1];

    ineqPolySys{idx}.noTerms = noTerms1;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineqSupp1;
    ineqPolySys{idx}.coef = coef1;
    idx = idx+1;
end

%2
noTerms2 = 3;
for k = 1:N
    ineqSupp2 = zeros(noTerms2, dimVar);
    ineqSupp2(1:2, 4:5) = eye(2);
    ineqSupp2(3, 5+k) = 1;
    coef2 = zeros(noTerms2, 1);
    coef2 = [-u(k); -u(k)^3; 1];

    ineqPolySys{idx}.noTerms = noTerms2;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 1;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineqSupp2;
    ineqPolySys{idx}.coef = coef2;
    idx = idx+1;
end


lbd = [-1; -1e10*ones(2, 1); 0.25; -1e10*ones(N+1, 1); -dEta*ones(N, 1)];
ubd = [1; 1e10*ones(2, 1); 0.25; 1e10*ones(N+1, 1); dEta*ones(N, 1)];

param.POPsolver = 'active-set';
param.relaxOrder = 2;
objPoly.degree = 1;
objPoly.noTerms = 1;
objPoly.dimVar = np+2*N;
objPoly.typeCone = 1;

pmin = zeros(np, 1);
pmax = zeros(np, 1);
for i = 1:np
    %min
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = 1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    %max
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = -1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end
    
[pmin, pmax]