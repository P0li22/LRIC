%% lab 8
clear all
close all
clc

% plant
Ts = 1/100;
s = tf('s');
z = tf('z', Ts);
Gcont = 4500/(s^2+16*s+4500);
Gz = c2d(Gcont, Ts, 'zoh');

% reference model
Mcont = (51.54*s+861.8)/(s^2+56.97*s+861.8);
Mz = c2d(Mcont, Ts, 'zoh');

% open loop experiment
N = 100;
r = unifrnd(-1, 2, N, 1);
y = lsim(Gz, r);
dEta = 0.01;
eta = unifrnd(-dEta, dEta, N, 1);
y_tilde = y+eta;

%s
Gm = zpk(minreal( Mz/(1-Mz) ));
s = lsim(Gm, r);

% opt problem
np = 7; dimVar = np+N; noTerms = 12;
a = 4;

ineqSupp = zeros(noTerms, dimVar);
A = zeros(noTerms, np);
B = zeros(noTerms, N);

% time invariant part of A
A(2:4, 1:3) = eye(3);
j = 1;
for i = 1:4
    A(4+j, 3+i) = 1;
    A(4+j+1, 3+i) = 1;
    j = j+2;
end

%time variant part of A and ineq definition
idx = 1;
for k = 4:N
    %fill B
    B = zeros(noTerms, N);
    j = 1;
    for i = 1:4
        B(5+j, k-i+1) = 1;
        j = j+2;
    end
    ineqSupp = [A, B];
    
    coef = zeros(noTerms, 1);
    for i = 1:4
        coef(i) = s(k-i+1);
    end
    j = 5;
    for i = 1:4
        coef(j) = -y_tilde(k-i+1);
        coef(j+1) = 1;
        j = j+2;
    end

    ineqPolySys{idx}.noTerms = noTerms;
    ineqPolySys{idx}.dimVar = dimVar;
    ineqPolySys{idx}.degree = 2;
    ineqPolySys{idx}.typeCone = -1;
    ineqPolySys{idx}.supports = ineqSupp;
    ineqPolySys{idx}.coef = coef;
    idx = idx+1;
end

lbd = [-1e10*ones(np, 1); -dEta*ones(N, 1)];
ubd = [1e10*ones(np, 1); dEta*ones(N, 1)];

pmin = zeros(np, 1);
pmax = zeros(np, 1);
for i = 1:np
    param.relaxOrder = 1;
    param.POPsolver = 'active-set';

    %min
    obj_supp = zeros(1, dimVar);
    obj_supp(i) = 1;
    obj_coef = 1;
    objPoly.noTerms = 1;
    objPoly.dimVar = dimVar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;

    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmin(i) = POP.xVectL(i);

    %max
    obj_supp = zeros(1, dimVar);
    obj_supp(i) = 1;
    obj_coef = -1;
    objPoly.noTerms = 1;
    objPoly.dimVar = dimVar;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = obj_supp;
    objPoly.coef = obj_coef;

    [~, ~, POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    pmax(i) = POP.xVectL(i);
end

[pmin, pmax]

%% central estimate, controller, step response
p = (pmax+pmin)/2;
Kc = tf([p(4) p(5) p(6) p(7)], [1 p(1) p(2) p(3)], Ts)
T = (Kc*Gz)/(1+Kc*Gz);
step(T);
hold on
step(Mcont);