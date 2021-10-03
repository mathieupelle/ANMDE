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

%% d) e) Computaation of the Fourier differentiation matrix

N_lst = 2:2:50; %Create array of elements to be analysed

error = zeros(length(N_lst),1); %Initialise error

%Inputs
v = @(x) exp(sin(pi*x)); %Test function we want to try
NST = false; %Whether or not to use negative sum trick
figure('Name','(e) - Differentiation matrix result for increasing N')

%Derive for each N
for n=1:length(N_lst)
    %Create the uniform spacing REMOVING THE ENDPOINT
    N = N_lst(n);
    N_step = 2/N;
    x = 0:N_step:2-N_step; %Notice the 2-N_step as endpoint
    D = get_FourierDiffMatrix(N,NST);
    %Compute the derivative once D has been calculated
    dvdx_fourier = D*v(x)'*pi;
    dvdx_analytic = pi*cos(pi*x).*exp(sin(pi*x));
    %Calculate analytical derivative and compare
    error(n,1) = sqrt(trapz(x,(dvdx_fourier - dvdx_analytic').^2));
    
    if N == 4 || N == 8 || N == 16 || N == 32
        plot(x, dvdx_fourier, 'DisplayName', sprintf('N = %d', N_lst(n)))
        hold on
        grid on
        xlabel('x')
        ylabel('$dv/dx$')
        legend('Location','best')
    end

end
plot(x, dvdx_analytic, '--k', 'DisplayName', 'Analytic')
saveas(gcf,'figures/e_derivative.eps','epsc')
% Do the same with NST
NST = true;
for n=1:length(N_lst)
    %Create the uniform spacing REMOVING THE ENDPOINT
    N = N_lst(n);
    N_step = 2/N;
    x = 0:N_step:2-N_step; %Notice the 2-N_step as endpoint
    D = get_FourierDiffMatrix(N,NST);
    %Compute the derivative once D has been calculated
    dvdx_fourier = D*v(x)'*pi;
    dvdx_analytic = pi*cos(pi*x).*exp(sin(pi*x));
    %Calculate analytical derivative and compare
    error(n,2) = sqrt(trapz(x,(dvdx_fourier - dvdx_analytic').^2));
end



%Get convergence error
p = polyfit(log(N_lst), log(error(:,1)), 1);
order = p(1);

figure('Name','(e)-Convergence using the fourier differenciation matrix')
loglog(N_lst, error(:,1), '-o'); hold on
loglog(N_lst, error(:,2), '-x')
legend('No NST','NST applied')
grid on; box on;
xlabel('N')
ylabel('$||dv/dx - D v_N||_{L^2}$')
saveas(gcf,'figures/e_convergence.eps','epsc')

figure('Name','(e) - NST results')
loglog(N_lst,abs(error(:,1)-error(:,2)))
grid on; box on;
xlabel('N')
ylabel('Error comparison')
saveas(gcf,'figures/e_NST.eps','epsc')

%% f) 


%Define all the functions and its derivatives analytically
syms x 
w0_l = -cos(pi*x);
w0_r = cos(pi*x);
w1_l = int(w0_l, x);
w1_r = int(w0_r, x);
w2_l = int(w1_l, x)-1/pi^2;
w2_r = int(w1_r, x)+1/pi^2;
w3_l = int(w2_l, x);
w3_r = int(w2_r, x);
w = [w0_l, w0_r; w1_l, w1_r; w2_l, w2_r; w3_l, w3_r];


NST = true;
N_lst = [2^2,2^3,2^4,2^5,2^6,2^7];
error = zeros(length(w)-1,length(N_lst));
figure('Name','f - Differentiation results with smoothness')
for n=1:length(N_lst)
    N = N_lst(n);
    N_step = 4/N;
    x = -2:N_step:2-N_step;
    %For each function, compute the derivative using D
    for f=1:length(w)-1
        D = get_FourierDiffMatrix(N,NST);
        
        %Divide between right and left part
        wl_l = matlabFunction(w(f,1));
        wl_r = matlabFunction(w(f,2));
        wu_l = matlabFunction(w(f+1,1));
        wu_r = matlabFunction(w(f+1,2));
        u = [wu_l(x(x<0)), wu_r(x(x>=0))];
        dwdx_fourier = D*u'*pi/2;
        dwdx_analytic = [wl_l(x(x<0)), wl_r(x(x>=0))];
        %error(f,n) = norm(dwdx_fourier' - dwdx_analytic);
        error(f,n) = sqrt(trapz(x,(dwdx_fourier' - dwdx_analytic).^2));
        %error(f,n) = L2norm(abs(dwdx_fourier - dwdx_analytic'), x);
        if n == length(N_lst)
            subplot(3,1,f)
            plot(x, dwdx_analytic, 'DisplayName', sprintf('$w_{%d}$', f-1));
            hold on; grid on; box on;
            plot(x, dwdx_fourier, 'DisplayName', sprintf('$v_{%d}$', f-1));
            legend('Location','best')
            ylabel('$w_i$')
            if f == 3
                ylim([-0.2,0.5])
            end
        end
    end
end
plot(x, u, 'DisplayName', sprintf('$w_{%d}$', f));
grid on
xlabel('x')
saveas(gcf,'figures/f_Derivatives.eps','epsc')


%Find order of convergence
ord = zeros(3,2);
for i=1:3
    ord(i,:) = polyfit(log10(N_lst),log10(error(i,:)),1);
end

figure('Name','f-Convergence of the different functions')
loglog(N_lst, error(1,:)/error(1,1),'-x' ,'DisplayName', [sprintf('$w_{%d}$', 0) + ": $O($"+num2str(round(ord(1,1),1))+"$)$"]);
hold on
loglog(N_lst, error(2,:)/error(2,1), '-x','DisplayName', [sprintf('$w_{%d}$', 1)+ ": $O($"+num2str(round(ord(2,1),1))+"$)$"]);
loglog(N_lst, error(3,:)/error(3,1),'-x' ,'DisplayName', [sprintf('$w_{%d}$', 2)+ ": $O($"+num2str(round(ord(3,1),1))+"$)$"]);
grid on
xlabel('N')
ylabel('$||w_i - D w_{i+1,N}||_{L^2}$')
legend('Location','best')
saveas(gcf,'figures/f_Convergence.eps','epsc')


%% g)

N_lst = [4,8,16,32,64,100,260,500,1000,2600,5000,10000];
%N_lst = 4:2:200;
%N_lst = [10]
error = zeros(length(N_lst),2);
time = zeros(length(N_lst),2);
v = @(x) exp(sin(pi*x));
NST = true;
% fftw('dwisdom',[]);
% fftw('planner','measure');
% X = rand(500,1);
% tic; fft(X); toc;
% fftinfo = fftw('dwisdom');
% fftw('dwisdom',fftinfo);
for n=1:length(N_lst)
    N = N_lst(n);
    N_step = 2/N;
    x = 0:N_step:2-N_step;
    
   
    
    
    k = [0:round(N/2-1), - round(N/2):-1];
    tic;
    dvdx_fft = ifft(1i*k.*fft(v(x)),'symmetric')*pi;
    time(n,2) = toc;
    
     
    D = get_FourierDiffMatrix(N,NST);
    
    tic;
    dvdx_fourier = D*v(x)'*pi;
    time(n,1) = toc; 


    dvdx_analytic = pi*cos(pi*x).*exp(sin(pi*x));
    error(n,1) = sqrt(trapz(x,(dvdx_fourier - dvdx_analytic').^2));
    %error(n,2) = sqrt(trapz(x,(real(dvdx_fft) - dvdx_analytic).^2));
    error(n,2) = sqrt(trapz(x,((dvdx_fft) - dvdx_analytic).^2));

end


figure('Name','g-Differentiation matrix vs FFT convergence')
loglog(N_lst, error(:,1),'-x', 'DisplayName', 'Differentiation matrix')
hold on
loglog(N_lst, error(:,2),'-x', 'DisplayName', 'FFT')
grid on
xlabel('N')
ylabel('$||e||_{L^2}$')
legend
saveas(gcf,'figures/g_convergence.eps','epsc')


%N vector for plotting convergence
N_con = 1000:50:10000;

figure('Name','g-Differentiation matrix vs FFT cpu time')
loglog(N_lst, time(:,1)*1000,'-x', 'DisplayName', 'Differentiation matrix')
hold on
loglog(N_lst, time(:,2)*1000,'-x', 'DisplayName', 'FFT')
loglog(N_con,(N_con.^2)*time(1,1)*1000/(N_con(1)^2),'--','Color','k','DisplayName','$N^2$')
loglog(N_con,N.*log10(N_con)*time(1,2)^1.4*1000/(N_con(1)*log10(N_con(1))),'-.','Color','k','DisplayName','$N \log N$')
grid on
xlabel('N')
%ylabel('CPU time [ms]')
legend('Location','best')
saveas(gcf,'figures/g_time.eps','epsc')
