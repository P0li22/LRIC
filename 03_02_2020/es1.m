clear all
close all
clc

load("data_exam_1A.mat");
u = input_data;
y = output_data;
np = 5;
N = 50;
dEta = 1;

%sparsePOP problem definition
noTerms = 9; dimVar = np+N; 

idx = 1;
for k=3:N
    ineqSupp = zeros(noTerms, dimVar);

    %fill time invariant part
    ineqSupp(2:6, 1:5) = eye(5);
    ineqSupp(8:9, 1:2) = eye(2);

    %fill time variant part
    ineqSupp(7:9, np+k-2:np+k) = flip(eye(3));

    coef = [y(k:-1:k-2)'; -u(k:-1:k-2)'; -1; -1; -1];

    ineqPolySys{idx}.supports = ineqSupp;
    ineqPolySys{idx}.coef = coef;
    ineqPolySys{idx}.noTerms = noTerms;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.typeCone = -1;

    idx = idx+1;
end

lbd = [-1e10*ones(np, 1); -dEta*ones(N, 1)];
ubp = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

pmin = zeros(np, 1);
pmax = zeros(np, 1);

objPoly.noTerms = 1;
objPoly.dimVar = dimVar;
objPoly.degree = 1;
objPoly.typeCone = 1;
for i = 1:np
    %min
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = 1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubp, param);
    pmin(i) = POP.xVectL(i);

    %max
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    objPoly.coef = -1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubp, param);
    pmax(i) = POP.xVectL(i);
end

PUI = [pmin, pmax]
pc = (pmin+pmax)/2