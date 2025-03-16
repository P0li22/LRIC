clear all
close all
clc

n = 2;
dEta = 0.1;
N = 100;
z = tf('z', 1);

load("data_exam_2A.mat");
r = Input;
y = Output;
np = 2*n+1;
M = (0.7009*z-0.5005)/(z^2-1.123*z+0.3234);
s = lsim(minreal(M/(1-M)), r);

% sparsePOP problem definition
noTerms = 3*n+3; dimVar = 2*n+N+1;

idx = 1;
for k=n+1:N
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

lbd = [-1e10*ones(2*n+1,1); -dEta*ones(N, 1)];
ubd = -lbd;
param.relaxOrder = 1;
param.POPsolver = 'active-set';

objPoly.typeCone = 1;
objPoly.dimVar = dimVar;
objPoly.degree = 1;
objPoly.noTerms = 1;

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

%% controller
p = (pmin+pmax)/2;
K = tf([p(3) p(4) p(5)], [1 p(1) p(2)], 1)