
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

%% Part h)

N = 6;
X = -1:0.01:1;

[Leg,Che,Tn] = deal(zeros(N,length(X)));
for x=1:length(X)
    Leg(:,x) = JacobiP(X(x),0,0,N)';
    Che(:,x) = JacobiP(X(x),-0.5,-0.5,N)';
end

%Scale cheby polys and create the legend
leg = [""];
for i=1:N
    Tn(i,:) = cheby_scale(Che(i,:),i-1);
    leg(i) = ["$P^{(0,0)} (N= $" + num2str(i-1) + ")"];
    leg_ch(i) = ["N=" + num2str(i-1)];
end

%Plotting
figure('Name','Legendre polys')
plot(X,Leg)
xlabel('x')
legend(leg)
grid on; box on;

figure('Name','Cheby polys')
plot(X,Tn)
xlabel('x')
legend(leg_ch,'Location','bestoutside'	)
grid on; box on;

P';
% L =[1, x, (3*x^2)/2 - 1/2, 0.5*(5*x^3-3*x), 1/8*(35*x^4-30*x^2+3)]
% L = [legendreP(0,x), legendreP(1,x), legendreP(2,x), legendreP(3,x), legendreP(4,x)]
C = [chebyshevT(0,x), chebyshevT(1,x), chebyshevT(2,x), chebyshevT(3,x), chebyshevT(4,x)]


%% F U N C T I O N S

function [P] = JacobiP(x, alpha, beta, N)

a = alpha; b = beta;
P = zeros(N,1);
P(1,1) = 1;
P(2,1) = 0.5*(a-b+x*(a+b+2));

for i=2:N-1
    n = i-1;
    a_l = 2*(n+a)*(n+b)/((2*n+a+b+1)*(2*n+a+b));
    a_c = (a^2-b^2)/((2*n+a+b+2)*(2*n+a+b));
    a_r = 2*(n+1)*(n+a+b+1)/((2*n+a+b+2)*(2*n+a+b+1));
    P(i+1,1) = ((a_c + x)*P(i) - a_l*P(i-1))/a_r;
end
end

function Tn = cheby_scale(pn,n)
    Tn = gamma(n+1)*gamma(0.5)/gamma(n+0.5) * pn;
end



