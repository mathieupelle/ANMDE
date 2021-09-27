clear all
close all
clc

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','docked')

%%

n = 3;
a = 0;
b = 0;


syms x
P = zeros(n,1);

P0 = 1;
P1 = 0.5*(a-b+x*(a+b+2));
P2 = (3*x^2)/2 - 1/2;
P3 =(x*(5*x^2 - 3))/2;
 

a_l = 2*(n+a)*(n+b)/((2*n+a+b+1)*(2*n+a+b));
a_c = a^2-b^2/((2*n+a+b+2)*(2*n+a+b));
a_r = 2*(n+1)*(n+a+b+1)/((2*n+a+b+2)*(2*n+a+b+1));

P4 = simplify(((a_c + x)*P3 - a_l*P2)/a_r)

%%
JacobiP(2, 0, 0, 1)

%% F U N C T I O N S


clear all
close all
clc

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','docked')

%%

N = 5;

a = 0;
b = 0;

x = 2;
P = zeros(N,1);
P(1,1) = 1;
P(2,1) = 0.5*(a-b+x*(a+b+2));

for i=2:N-1
    n = i-1;
    a_l = 2*(n+a)*(n+b)/((2*n+a+b+1)*(2*n+a+b));
    a_c = a^2-b^2/((2*n+a+b+2)*(2*n+a+b));
    a_r = 2*(n+1)*(n+a+b+1)/((2*n+a+b+2)*(2*n+a+b+1));
    P(i+1,1) = ((a_c + x)*P(i) - a_l*P(i-1))/a_r;
end

P'
L =[1, x, (3*x^2)/2 - 1/2, 0.5*(5*x^3-3*x), 1/8*(35*x^4-30*x^2+3)]

%%
JacobiP(2, 0, 0, 1)

%% F U N C T I O N S

function [P] = JacobiP(x, alpha, beta, n)    
    P = jacobiP(n, alpha, beta, x);
end



