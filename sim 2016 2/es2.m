clear all
close all
clc

load("data_exam_2A_dg.mat");
u = Input;
y = Output;
dEta = 0.07; 
dEps = 0.07;
np = 5;
N = 30;

%sparsePOP problem definition
%1
noTerms1 = 8; dimVar = 5+2*N; 

idx = 1;
for k=3:N
    ineqSupp1 = zeros(noTerms1, dimVar);
    ineqSupp1(2:4, 1:3) = eye(3);
    ineqSupp1(6:8, 1:3) = eye(3);
    ineqSupp1(5:6, np+k-1:np+k) = flip(eye(2));
    ineqSupp1(7:8, np+N+k-1:np+N+k) = flip(eye(2));
    coef1 = [y(k:-1:k-1); -u(k:-1:k-1); -1; -1; 1; 1];

    ineqPolySys{idx}.noTerms = noTerms1;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineqSupp1;
    ineqPolySys{idx}.coef = coef1;

    idx = idx+1;
end

%2
noTerms2 = 4; 
ineqSupp2 = zeros(noTerms2, dimVar);
ineqSupp2(1:2, 2:3) = eye(2);
ineqSupp2(3, 1) = 1;
ineqSupp2(3:4, 4) = [1; 1];
coef2 = [1, 1, -1, -1]';

ineqPolySys{idx}.noTerms = noTerms2;
ineqPolySys{idx}.dimVar = dimVar;
ineqPolySys{idx}.degree = 2;
ineqPolySys{idx}.typeCone = -1;
ineqPolySys{idx}.supports = ineqSupp2;
ineqPolySys{idx}.coef = coef2;

idx = idx+1;

%3
noTerms3 = 2;
ineqSupp3 = zeros(noTerms3, dimVar);
ineqSupp3(1, 3) = 1;
ineqSupp3(2, 2) = 1;
ineqSupp3(2, 5) = 1;
coef3 = [1; 1];

ineqPolySys{idx}.noTerms = noTerms3;
ineqPolySys{idx}.dimVar = dimVar;
ineqPolySys{idx}.degree = 2;
ineqPolySys{idx}.typeCone = -1;
ineqPolySys{idx}.supports = ineqSupp3;
ineqPolySys{idx}.coef = coef3;

lbd = [-1e10*ones(3, 1); 5; 0.1; -dEta*ones(N, 1); -dEps*ones(N, 1)];
ubd = [1e10*ones(3, 1); 9; 0.5; dEta*ones(N, 1); dEps*ones(N, 1)];

param.relaxOrder = 1;
param.POPsolver = 'active-set';

objPoly.degree = 1;
objPoly.typeCone = 1;
objPoly.noTerms = 1;
objPoly.dimVar = dimVar;
pmin = zeros(3, 1);
pmax = zeros(3, 1);
for i = 1:3
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

%% sim
G = tf([pc(2:3)'],[1 pc(1)],0.1);

y_sim = lsim(G,u);

figure
plot(y,'b')
hold on 
plot(y_sim,'r--')
grid on