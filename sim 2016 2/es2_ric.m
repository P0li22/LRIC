clear all
close all
clc

load data_exam_2A_dg.mat


utilde = Input;
ytilde = Output;

Deta = 0.07;
Deps= 0.07;

N = 30;
n = 1;


dimvar = 3+2*N+2;

c = 1;


for k = n+1:N

    support = zeros(8,dimvar);
    support(2,3+k) = 1;
    support(3,1) = 1;
    support(4,1) = 1;
    support(4,3+k-1) = 1;
    support(5,2) = 1;
    support(6,2) = 1;
    support(6,3+N+k) = 1;
    support(7,3) = 1;
    support(8,3) = 1;
    support(8,3+N+k-1) = 1;


    ineqPolySys{c}.noTerms = 8;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.dimVar = dimvar;
    ineqPolySys{c}.typeCone = -1;

    ineqPolySys{c}.supports = support;
    ineqPolySys{c}.coef = [ytilde(k);-1;ytilde(k-1);-1;-utilde(k);1;-utilde(k-1);1];

    c = c+1;



end







ineqPolySys{c}.noTerms = 4;
ineqPolySys{c}.dimVar = dimvar;
ineqPolySys{c}.degree = 2;
ineqPolySys{c}.typeCone = -1;

support = zeros(4,dimvar);

support(1,1) = 1;
support(1,3+2*N+1) = 1;
support(2,2) = 1;
support(3,3) = 1;
support(4,3+2*N+1) = 1;


ineqPolySys{c}.supports = support;
ineqPolySys{c}.coef = [1;-1;-1;1];



c = c+1;





ineqPolySys{c}.noTerms = 2;
ineqPolySys{c}.dimVar = dimvar;
ineqPolySys{c}.degree = 2;
ineqPolySys{c}.typeCone = -1;

support = zeros(2,dimvar);

support(1,2) = 1;
support(1,3+2*N+2) = 1;
support(2,3) = 1;


ineqPolySys{c}.supports = support;
ineqPolySys{c}.coef = [1;1];






objPoly.noTerms = 1;
objPoly.dimVar = dimvar;
objPoly.degree = 1;
objPoly.typeCone = 1;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

ubd = [1e10*ones(3,1);Deta*ones(N,1);Deps*ones(N,1);9;0.5];
lbd = [-1e10*ones(3,1);-Deta*ones(N,1);-Deps*ones(N,1);5;0.1];



theta = zeros(3,2);



for p = 1:3


    supp = zeros(1,dimvar);
    supp(1,p) = 1;
    objPoly.supports = supp;

    objPoly.coef = 1;
    [a,b,POP] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
    theta(p,1) = POP.objValueL;



    objPoly.coef = -1;
    [a,b,POP] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
    theta(p,2) = -POP.objValueL;



end

theta


theta_c = (theta(:,1)+theta(:,2))/2;
theta_c = theta_c';



G = tf(theta_c(2:3),[1 theta_c(1)],0.1);

y_sim = lsim(G,utilde);

figure
plot(ytilde,'b')
hold on 
plot(y_sim,'r--')
grid on









