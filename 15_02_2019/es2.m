clear all
close all
clc

z = tf('z', 1)

load("data_exam_2A.mat");
dEta = 0.1;
n = 2;
N = 100;
r = Input;
y = Output;
M = (0.7009*z-0.5005)/(z^2-1.123*z+0.3234);
s = lsim(minreal(M/(1-M)), r);

%sparsePOP problem definition
np = 2*n+1; noTerms = 3*n+3; dimVar = np+N;

idx = 1;
for k = n+1:N
    ineqSupp = zeros(noTerms, dimVar);
    ineqSupp(2:np+1, 1:np) = eye(np);
    ineqSupp(np+2:end, n+1:np) = eye(n+1);
    ineqSupp(np+2:end, np+k-n:np+k) = flip(eye(n+1));
    coef = [s(k:-1:k-n); -y(k:-1:k-n); ones(n+1, 1)];

    ineqPolySys{idx}.noTerms = noTerms;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.supports = ineqSupp;
    ineqPolySys{idx}.coef = coef;
    idx = idx+1;
end

lbd = [-1e10*ones(np, 1); -dEta*ones(N, 1)];
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

pmin = zeros(np, 1);
pmax = zeros(np, 1);

objPoly.noTerms = 1;
objPoly.dimVar = dimVar;
objPoly.typeCone = 1;
objPoly.degree = 1;
for i=1:np
    objSupp = zeros(1, dimVar);
    objSupp(i) = 1;
    objPoly.supports = objSupp;
    %min
    objPoly.coef = 1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    %max
    objPoly.coef = -1;
    [aa, bb, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end

[pmin, pmax]
pc = (pmin+pmax)/2
K = tf([pc(n+1:np)'], [pc(1:n)'], 1)