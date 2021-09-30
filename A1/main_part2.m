
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
leg_ch = [""];
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

% Function to interpolate
u = @(x) 1./(2-cos(pi*x)); % x on [0,2]

% Order --> number of coefficients
K = 200;

% N+1 number of points
N_lst = [10, 40, 80, 100, 200];
%N_lst = [2:1:40];
error = zeros(length(N_lst), 1); %error vector
%K = 10;
%N_lst = [8, 9, 10];

f1 = figure('Name' , 'Comparison');
f2 = figure('Name', 'Decay');
%Loop over number of points
for n=1:length(N_lst)
    N = N_lst(n);
    [x_basis, w] = JacobiGQ(0, 0, N); %Gauss-Legendre nodes and weights
    
    %Basis function
    phi = zeros(K, N+1); 
    for i=1:N+1
        phi(:,i) = JacobiP(x_basis(i), 0, 0, K);
    end

    % Scaling of nodes and exclude endpoints (Gauss-Type)
    a = 0;
    b = 2;
    x = (b-a)/2.*x_basis + (b+a)/2; %scaling
    
    % Evaluate coefficients using function at scaled node locations
    f = zeros(K,1);
    for k=1:K
        gamma = 0;
        for j=1:N+1
            gamma = gamma + phi(k, j)^2*w(j);
        end
        f(k,1) = 1/gamma*sum(u(x).*phi(k,:)'.*w); %coefficients
    end
    
    e = (u(x)'- f'*phi)'; %error
    error(n,1) = L2norm(e, x_basis); %L2 norm
    
    set(0, 'CurrentFigure', f1)
    plot(x, f'*phi, '--o', 'DisplayName', sprintf('N = %d', N))
    hold on
    
    % Decaying coefficients
    if N == N_lst(end)
        set(0, 'CurrentFigure', f2)
        plot(1:K, f)
        grid on
        ylim([0,max(f)]);
        xlabel('k')
        ylabel('$\bar{f}_k$')
        axes('position',[.54 .430 .35 .35])
        box on
        stem(1:20, f(1:20))
        %axis tight
        xlabel('k')
        ylabel('$\bar{f}_k$')
        grid on
    end
    
end
set(0, 'CurrentFigure', f1)
x = linspace(0,2,50);
plot(x, u(x), '--k', 'DisplayName', 'Analytical')
grid on
legend('Location', 'Best')
xlabel('x')
ylabel('u(x)')


figure('Name', 'Error')
loglog(N_lst, error, '-x')
grid on
xlabel('N')
ylabel('Error, $\| u(x)-\mathcal{P}_Nu(x)  \|_{L2}$ ')

%% part j)

N = 5;
x = JacobiGL(0, 0, N);

V = zeros(N+1, N+1);
for i=1:N+1
    V(:,i) = JacobiP(x(i), 0, 0, N+1);
end

x = linspace(-1, 1, 100);

h = zeros(N+1,length(x));
for i=1:length(x)
    phi = JacobiP(x(i), 0, 0, N+1);
    h(:,i) = (V)^(-1)*phi;
end

figure
for i=1:N
    plot(x, h(i,:), 'DisplayName', sprintf('$h_{%d}$', i-1))
    hold on
end
grid on
xlabel('x')
ylabel('$h_i$')
legend('Location','EastOutside')
x = JacobiGL(0, 0, N);
scatter(x, zeros(length(x)), 'ko', 'HandleVisibility', 'off')

%% part k)

% Function to differentiate
v = @(x) exp(sin(pi*x)); % x on [0,2]

N_lst = [2:5:100,100:100:1000];
%N_lst = [5,10,20,40];
error = zeros(length(N_lst), 1);
L2 = zeros(length(N_lst), 1);
figure;
for n=1:length(N_lst)
    N = N_lst(n);
    x_basis = JacobiGL(0, 0, N); % Get the nodes (Gauss Lobato), using both endpoints (not idiot)

    [V,Vx] = deal(zeros(N+1, N+1)); 
    % Loop through each node, and evaluate Vandermonde and its derivative (not
    % idiot, probably) (using Legendre polynomials)
    for i=1:N+1
        V(:,i) = JacobiP(x_basis(i), 0, 0, N+1);
        Vx(:,i) = gradJacobiP(x_basis(i), 0, 0, N+1);
    end

    %Create differenciation matrix
    D = Vx'*(V')^(-1);
    %We have evaluated the Vandermome matrices from -1 to +1, and now we need
    %to map those nodes to the domain of our function
    x = (x_basis + 1)';
    %Evaluate derivative
    dvdx = D*v(x)';
    dvdx_analytic = pi*cos(pi*x).*exp(sin(pi*x));
    error(n, 1) = L2norm((dvdx - dvdx_analytic'), x_basis);

    if n <= length(N_lst)
        plot(x,dvdx, '-x', 'DisplayName', sprintf('N = %d',N))
        hold on
        
    end
end
x = linspace(0,2,100);
dvdx_analytic = pi*cos(pi*x).*exp(sin(pi*x));
plot(x,dvdx_analytic, '--k', 'DisplayName', 'Analytical')
grid on
xlabel('x')
ylabel('$\frac{dv}{dx}$', 'FontSize', 18)
legend('Location', 'Best')

figure('Name', 'Error')
loglog(N_lst, error, '-x')
grid on
xlabel('N')
ylabel('Error, $\| \frac{dv(x)}{dx}-\mathcal{D}v(x)  \|_{L2}$ ')

%% part l)

syms x
u1 = @(x) 0*x+1;
u2 = @(x) x;
u3 = @(x) x.^2;
u4 = @(x) x.^3;
u5 = @(x) sin(x);

L2_ana1 = sqrt(double(simplify(int(0*x+1, 0, 2))));
L2_ana2 = sqrt(double(simplify(int(x^2, 0, 2))));
L2_ana3 = sqrt(double(simplify(int(x^4, 0, 2))));
L2_ana4 = sqrt(double(simplify(int(x^6, 0, 2))));
L2_ana5 = sqrt(double(simplify(int(sin(x)^2, 0, 2))));

N = 1:1:5;
L2_num = zeros(length(N),2);
for n=1:length(N)
    L2_num(n,1) = MatrixBasedInt(u1, 0, 2, N(n));
    L2_num(n,2) = MatrixBasedInt(u2, 0, 2, N(n));
    L2_num(n,3) = MatrixBasedInt(u3, 0, 2, N(n));
    L2_num(n,4) = MatrixBasedInt(u4, 0, 2, N(n));
    L2_num(n,5) = MatrixBasedInt(u5, 0, 2, N(n));
end

figure
plot(N, abs(L2_num(:,1)-L2_ana1), '-o')
hold on 
plot(N, abs(L2_num(:,2)-L2_ana2), '-x')
plot(N, abs(L2_num(:,3)-L2_ana3), '-x')
plot(N, abs(L2_num(:,4)-L2_ana4), '-x')
plot(N, abs(L2_num(:,5)-L2_ana5), '-x')
% plot(N, L2_ana1*ones(length(N)), '--', 'Color', [0 0.4470 0.7410])
% plot(N, L2_ana2*ones(length(N)), '--', 'Color', [0.8500 0.3250 0.0980])
grid on
xlabel('N')
ylabel('Error, $|\sqrt{\mathbf{f^T}\mathcal{M}\mathbf{f}}- \| u(x)  \|_{L2}|$')
legend('$u(x)=1$', '$u(x)=x$', '$u(x)=x^2$', '$u(x)=x^3$','$u(x)=sin(x)$')

% up to polynomial in basis
%% F U N C T I O N S

function [L2] = MatrixBasedInt(fun, a, b, N)

    x = JacobiGL(0, 0, N); % Get the nodes (Gauss Lobato), using both endpoints (not idiot)
    V = zeros(N+1, N+1); 

    for i=1:N+1
        P = JacobiP(x(i), 0, 0, N+1);
        V(:,i) = OrthoJacobiP(P, 0, 0);
    end
    size(V)
    M = (V'*V)^(-1); % Mass matrix
    x = (b-a)/2.*x + (b+a)/2; % scaling
    f = fun(x);
    L2 = sqrt(f'*M*f); 

end


