clear all
close all
clc

load("data_exam_1A.mat");
u = Input;
y = Output;

N = 50;
dEta = 1;

%sparsePOP problem formulation
noTerms = 8; dimVar = 4+N; np = 4;

idx = 1;
for k=3:N 
    ineqSupp = zeros(noTerms, dimVar);
    ineqSupp(2:5, 1:4) = eye(4);
    ineqSupp(7:8, 1:2) = eye(2);
    ineqSupp(6:8, np+k-2:np+k) = flip(eye(3));

    coef = [y(k:-1:k-2); -u(k-1:-1:k-2); -1*ones(3, 1)];

    ineqPolySys{idx}.noTerms = noTerms;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineqSupp;
    ineqPolySys{idx}.coef = coef;

    idx = idx+1;
end

lbd = [-1e10*ones(4, 1); -dEta*ones(N, 1)];
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

objPoly.noTerms = 1;
objPoly.dimVar = dimVar;
objPoly.degree = 1;
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
pc = (pmin+pmax)/2