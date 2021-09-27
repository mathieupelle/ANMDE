
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

%% Part i)

u = @(x) 1./(2-cos(pi*x)); % x on [0,2]

K = 200;

N_lst = [10, 40, 80, 100, 200];
error = zeros(length(N_lst), 1);

f1 = figure('Name' , 'Comparison');
f2 = figure('Name', 'Decay');
for n=1:length(N_lst)
    N = N_lst(n);
    [x, w] = JacobiGQ(0, 0, N);

    phi = zeros(K, N+1);
    for i=1:N+1
        phi(:,i) = JacobiP(x(i), 0, 0, K);
    end

    dx = 2/(N+2);
    x = 0+dx:dx:2-dx;
    x = x';

    f = zeros(K,1);
    for k=1:K
        gamma = 0;
        for j=1:N+1
            gamma = gamma + phi(k, j)^2*w(j);
        end
        f(k,1) = 1/gamma*sum(u(x).*phi(k,:)'.*w);
    end
    
    error(n,1) = mean(abs(u(x)'- f'*phi));
    
    set(0, 'CurrentFigure', f1)
    plot(x, f'*phi)
    hold on
    
    if N == N_lst(end)
        set(0, 'CurrentFigure', f2)
        plot(1:K, f)
        hold on
        grid on
        xlabel('k')
        ylabel('$\bar{f}_k$')
    end
    
end
set(0, 'CurrentFigure', f1)
plot(x, u(x), 'x')
grid on


figure('Name', 'Error')
loglog(N_lst, error)
grid on
xlabel('N')
ylabel('error')

%% part j)

N = 6;
x = JacobiGL(0, 0, N);

V = zeros(N+1, N+1);
for i=1:N+1
    V(:,i) = JacobiP(x(i), 0, 0, N+1);
end

%x = linspace(-1, 1, 100);

h = zeros(N+1,length(x));
for i=1:length(x)
    phi = JacobiP(x(i), 0, 0, N+1);
    h(:,i) = (V')^(-1)*phi;
end

figure
for i=1:N
    plot(x, h(i,:), 'DisplayName', sprintf('$h_{%d}$', i-1))
    hold on
end
grid on
xlabel('x')
ylabel('$h_i$')
%legend
x = JacobiGL(0, 0, N);
scatter(x, zeros(length(x)))


%% part k)

% N = 6;
% x = JacobiGL(0, 0, N);
% 
% V = zeros(N+1, N+1);
% for i=1:N+1
%     V(:,i) = JacobiP(x(i), 0, 0, N+1);
% end

gradJacobiP(1, 0, 0, 3)



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

function [dPdx] = gradJacobiP(x, alpha, beta, N)
    
    dPdx = zeros(N, 1);
    P = JacobiP(x, alpha+1, beta+1, N);
    for i=1:N
        n = i-1;
        if i == 1
            dPdx(i,1) = 0;
        else
            dPdx(i,1) = sqrt(n*(n+alpha+beta+1))*P(i-1,1);
        end
    end
    
end



