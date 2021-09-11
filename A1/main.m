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

%% a)

x_min = 0;
x_max = 2;

syms x n
f = 1/(2-cos(pi*x));
cn = 1/2*int(cos(n*pi*x)/(2-cos(pi*x)), x, x_min, x_max);

N_max = 40;
N_lst = -N_max:N_max;
cn_val = zeros(length(N_lst), 1);
for j = 1:length(N_lst)
    %cn_val(j, 1) = double(subs(cn, n, N_lst(j)));
    cn_val(j,1) = (2-sqrt(3))^abs(N_lst(j))/sqrt(3);
end

N_lst_pos = N_lst(floor(N_max+1):end);
cn_val_pos = cn_val(floor(N_max+1):end);
p = polyfit(N_lst_pos, log(cn_val_pos), 1);
decay = p(1);

figure
semilogy(N_lst_pos, cn_val_pos,'k')
hold on
semilogy(N_lst_pos, exp(p(2))*exp(p(1)*N_lst_pos))
grid on
xlabel('n')
ylabel('$c_n$')


x_val = linspace(x_min, x_max, 100);

figure;
y_analytic = double(subs(f, x, x_val))';
plot(x_val, y_analytic, '--k')

N = 1:35;
y = zeros(length(x_val), length(N));
error = zeros(length(N),1);
for n=1:length(N)
    N_lst = -N(n):N(n);
    coef = cn_val(N_max-n-1:N_max+n-1);
    for j=1:length(x_val)
        yt = 0;
        for k=1:length(N_lst)
            yt = yt + coef(k)*exp(1i*N_lst(k)*x_val(j)*pi); 
        end
        y(j,n) = yt;   
    end
    error(n,1) = norm(y_analytic - abs(y(:,n)));
    plot(x_val, abs(y(:,n)), 'DisplayName', sprintf('N = %d', N(n)));
    hold on
end
xlabel('x')
ylabel('y')
grid on
%legend;

figure;
loglog(N, error,'-o');
grid 
xlabel('N')
ylabel('error')
grid on


%% b)

N_lst = 2.^(2:1:6);
f = @(x) 1./(2-cos(pi*x));

figure;
for j=1:length(N_lst)
    N = N_lst(j);
    x = linspace(x_min, x_max, N);
    u = f(x);
    cn = fft(u)/N;
    cn = fftshift(cn);
    n = -N/2:1:N/2-1;
    Pn = zeros(N,1);
    for g=1:N
        pt = 0;
        for k=1:length(n)
            pt = pt + cn(k)*exp(1i*n(k)*x(g)*pi); 
        end
        Pn(g,1) = pt;   
    end
    plot(x, abs(Pn), 'DisplayName', sprintf('N = %d', N));
    hold on;
end

x_analytic = linspace(x_min, x_max, 100);
plot(x_analytic, f(x_analytic), '--k', 'DisplayName', 'Analytical');
grid on;
legend;

%% c)

N = 6;
j_lst = 0:N-1;
x = linspace(0, 2*pi, 50);
h = zeros(N,length(x));

figure;
for j=1:length(j_lst)
    xj = 2*pi/N*j_lst(j);
    for i=1:length(x)
        h(j,i) = 1/N*sin(N/2*(x(i)-xj))*cot(0.5*(x(i)-xj));
    end
    plot(x, h(j,:), 'DisplayName', sprintf('$h_%d$', j));
    hold on
end
grid on
xlabel('x')
ylabel('$h_j$')
legend
