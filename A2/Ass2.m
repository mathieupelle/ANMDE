%% Advanced Numerical Methods for Differential Equations - Assignment 2
clear variables; close all; clc; beep off;

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
%set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','docked')
%% Question (a) - Legendre Tau Method (LTM)

%Input: 
N = [2 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9]; %how many basis functions?
eps = [0.1 0.01];% 0.001]; %Function parameter
error = zeros(length(N),length(eps));
X = 0:0.005:1;


%Call the solver for each eps value and number of basis
for i=1:length(eps)
    u_ana = (exp(-X/eps(i))+X-1-X.*exp(-1/eps(i)))./(exp(-1/eps(i))-1); %Get analytical solution
    for j=1:length(N)
        ux = LTM_a(N(j),eps(i)); %Get numerical solution
        
        %Compute error with L2 norm
        error(j,i) = sqrt(trapz(X,(u_ana-ux').^2));     
       
    end
end

figure('Name','a-Legendre Tau Method solution')
plot(X,ux)
hold on
plot(X,u_ana)
legend('Legendre Tau Method','Analytical')
xlim([0 1])
xlabel('$x$')
ylabel('$u$')

figure('Name','a-Legende Tau Method convergence')
loglog(N,error,'-x')
grid on; box on;
legend('$\epsilon = 0.1$','$\epsilon = 0.01$','$\epsilon = 0.001$')
xlabel('Number of basis $N$')
ylabel('Error L2 norm')

%% Question (a2) - Legendre Collocation Method


N=2^4;
eps = 0.01;
X = 0:0.005:1;
u_ana = (exp(-X/eps)+X-1-X.*exp(-1/eps))./((exp(-1/eps))-1); %Get analytical solution

for i=1:length(N)
    Ni = N(i);
    x_basis = JacobiGL(0, 0, Ni); % Get the nodes (Gauss Lobato), using both endpoints (not idiot)

    [V,Vx] = deal(zeros(Ni+1, Ni+1)); 
    % Loop through each node, and evaluate Vandermonde and its derivative (not
    % idiot, probably) (using Legendre polynomials)
    for i=1:Ni+1
        V(:,i) = JacobiP(x_basis(i), 0, 0, Ni+1);
        Vx(:,i) = gradJacobiP(x_basis(i), 0, 0, Ni+1);
    end

    %Create differenciation matrix
    D = Vx'*(V')^(-1);
    
    %Build differentation operand
    L = 4*eps*D^2+2*D;
    
    %Create RHS vector 
    f = -1*ones(Ni+1,1);
    
    %Impose BCs
    %L(1,:) = [1,zeros(1,Ni)];
    L(end-1,:) = [1,zeros(1,Ni)];
    L(end,:) = [zeros(1,Ni),1];
    %f(1) = 0;
    f(end-1) = 0;
    f(end) = 0;
    
    %Solve the system
    u = L\f;
end

figure
plot((x_basis+1)/2,u);
hold on
plot(X,u_ana)

%%  Question (b) - Irrotational flow around a cylinder 

%Number of points
Ni = 35;  %odd for even N+1

% Create coordinates array
r1 = 1; 
a = 3;
r2 = a*r1;

N_step = 2*pi/(Ni+1);
theta = 0:N_step:2*pi-N_step;

%Create differentation matrices
% For radial:
x_basis = JacobiGL(0, 0, Ni); % Get the nodes (Gauss Lobato), using both endpoints (not idiot)
r = x_basis*(r2-r1)/2+(r2+r1)/2;

%Create meshgrid
[R,TH] = meshgrid(r,theta);

%Vectorize mesghrid
R = R(:); TH = TH(:);

[V,Vx] = deal(zeros(Ni+1, Ni+1));
% Loop through each node, and evaluate Vandermonde and its derivative (not
% idiot, probably) (using Legendre polynomials)
for i=1:Ni+1
    V(:,i) = JacobiP(x_basis(i), 0, 0, Ni+1);
    Vx(:,i) = gradJacobiP(x_basis(i), 0, 0, Ni+1);
end

%Create differenciation matrix
Dr = Vx'*(V')^(-1);
Dr = Dr*(2/(r1*(a-1)))^-1;

%Get differentation matrix with Fourier basis
Dth = get_FourierDiffMatrix(Ni+1,false);
Dthth = Dth*Dth;

%Create them in 2D
I = eye(Ni+1);
DR = kron(Dr,I); 
DTT = kron(I,Dthth);

%Construc differentation operand
L = inv(diag(R)) * DR * (diag(R)*DR) + inv(diag(R)^2) * DTT;

%Build RHS
f = zeros((Ni+1)*(Ni+1),1);

%Plug boundary conditions
V_inf = 1;
f(1:Ni+1) = 2*V_inf*r1*cos(theta);
f(end-(Ni):end) = V_inf*(r2+r1^2/r2)*cos(theta); 
%Modify differentiation operand
L(1:Ni+1,1:Ni+1) = eye(Ni+1); L(1:Ni+1,Ni+2:end) = 0;
L(end-Ni:end,end-Ni:end) = eye(Ni+1); L(end-Ni:end,1:end-Ni-1) = 0;

%Solve the system
u = L\f;
u_reshaped = reshape(u(:),Ni+1,Ni+1);

%Compute the error
[R,TH] = meshgrid(r,[theta, 2*pi]);

u_reshaped = vertcat(u_reshaped, u_reshaped(1,1:Ni+1));
u_th = V_inf*(R+r1^2*R.^-1).*cos(TH);
err = abs((u_reshaped-u_th));

error = (u_reshaped-u_th).^2;
L2 = sqrt(trapz(theta,trapz(r,error(1:end-1,:))));

%Plotting
figure('Name','Error in logscale')
%contourf(R.*cos(TH),R.*sin(TH),reshape(u(:),Ni+1,Ni+1));
contourf(R.*cos(TH),R.*sin(TH),log10(err));
colorbar

figure('Name','Velocity potential')
subplot(121)
contourf(R.*cos(TH),R.*sin(TH),u_reshaped);
colorbar
subplot(122)
contourf(R.*cos(TH),R.*sin(TH),u_th);
colorbar

%%  Question (c) - Solving the KdV equation

saving_hist = 1;
conservation = 1;
dealiasing = 0;

Ni_lst = 2.^(3:7) + 1;
Ni_lst = 111;
finalL2norm = zeros(length(Ni_lst),1);
for n=1:length(Ni_lst)
    Ni = Ni_lst(n);
    txt = ['=> Computing for N = ',num2str(Ni)]; disp(txt)
    c = 1;
    x0 = 0;
    x_limits = [-2*pi, 6*pi];
    end_time = 10;

    [u_hist, u_ana, errors, x, time, ~] = RK4_KdV(Ni, c, x0, x_limits, end_time, saving_hist, conservation, dealiasing);
    
    finalL2norm(n, 1) = sqrt(abs(trapz((u_hist(:,end) - u_ana(:,end)).^2, x)));

    if n==length(Ni_lst)
        if saving_hist ==1

%             [X,T] = meshgrid(x,time);
%             figure('Name', 'Spectral method')
%             surf(X,T,u_hist')
%             xlabel('x')
%             ylabel('t')
%             zlabel('$\overline{u}(x,t)$')
% 
%             figure('Name', 'Soliton')
%             surf(X,T,u_ana')
%             xlabel('x')
%             ylabel('t')
%             zlabel('u(x,t)')
% 
%             figure('Name', 'Absolute error')
%             surf(X,T,abs(u_ana'-u_hist'))
%             xlabel('x')
%             ylabel('t')
%             zlabel('$|\overline{u}(x,t) - u(x,t)|$')

            figure('Name', 'L2 norm in time')
            plot(time, errors.L2norm)
            xlabel('Time')
            ylabel('$||u-\mathcal{I}_Nu||_{L2}$')
            grid on
        end
        figure('Name', 'Final timestep')
        plot(x, u_ana(:,end), 'DisplayName', 'Analytic')
        hold on
        plot(x, u_hist(:,end), '--.k', 'DisplayName', 'spectral method')
        xlabel('x')
        ylabel('u')
        legend
        grid on
    end
end


% figure('Name', 'Convergence')
% plot(x, finalL2norm)
% grid on
% xlabel('N')
% ylabel('$||u-\mathcal{I}_Nu||_{L2}$')

% if saving_hist == 1
    figure
    hold on
    for i = 1:1000:length(time)
        plot(x,abs(u_ana(:,i)-u_hist(:,i)), '-r')
        hold on
        %plot(x,u_ana(:,i), '--ok')
        %ylim([0 max(u_ana(:,1))*1.02]);
        pause(0.005)
        clf
    end
% end

%% Question (d) - conservation of quantities

saving_hist = 1;
conservation = 1;
dealiasing = 0;

c_lst = [0.25, 0.5, 1];
leg = {'c = 0.25', 'c = 0.5', 'c = 1.0'};
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
for i=1:length(c_lst)
    
    Ni = 151;
    x0 = 0;
    c = c_lst(i);
    x_limits = [-2*pi, 4*pi];
    end_time = 10;
    
    [u_hist, ~, errors, x, time, quant] = RK4_KdV(Ni, c, x0, x_limits, end_time, saving_hist, conservation, dealiasing);
    
    set(0, 'CurrentFigure', f1)
    semilogy(time, errors.norm2)
    grid on
    hold on
    xlabel('Time')
    ylabel('$||u-\mathcal{I}_Nu||_2$')
    legend(leg, 'Location', 'Best')
    
    set(0, 'CurrentFigure', f2)
    semilogy(time, errors.normInf)
    grid on
    hold on
    xlabel('Time')
    ylabel('$||u-\mathcal{I}_Nu||_\infty$')
    legend(leg, 'Location', 'Best')

    set(0, 'CurrentFigure', f3)
    err = abs(quant.M - quant.M_ana);
    semilogy(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Error in Mass')
    legend(leg, 'Location', 'Best')
    
    set(0, 'CurrentFigure', f4)
    err = abs(quant.V - quant.V_ana);
    semilogy(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Error in Momentum')
    legend(leg, 'Location', 'Best')
    
    
    set(0, 'CurrentFigure', f5)
    err = abs(quant.E - quant.E_ana);
    semilogy(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Error in Energy')
    legend(leg, 'Location', 'Best')
    
end
    figure
    hold on
    for i = 1:1000:length(time)
        plot(x,u_hist(:,i), '-r')
        hold on
        ylim([-0.5 3]);
        pause(0.005)
        clf
    end
%% Question (d) - error for different velocities

x0 = 0;
c_lst = [0.25, 0.5, 1];
x = linspace(0, 4*pi,10000);
soli = 0.5*c*sech(0.5.*sqrt(c_lst).*(x'-x0)).^2;
leg = {'c = 0.25', 'c = 0.5', 'c = 1.0', 'c = 5.0'};

figure('Name','Soliton  for domain size study')
semilogy(x,soli)
legend(leg, 'Location', 'Best')
xlabel('Position $x$')
grid on
ylabel('Soliton $u(x,t=0)$')

saving_hist = 1;
conservation = 1;
dealiasing = 0;
Ni_lst = 2.^(2:8) + 1;

x_limits = [-4*pi, 4*pi];
end_time = 0.1;

finalL2norm = zeros(length(Ni_lst),length(c_lst));
[err2,errinf,errl2] = deal(zeros(length(Ni_lst),length(c_lst)));
for c=1:length(c_lst)
    
    parfor n=1:length(Ni_lst)
        Ni = Ni_lst(n);
        txt = ['=> Computing for N = ',num2str(Ni)];
        disp(txt)

        
        [u_hist, u_ana, errors, x, ~, ~] = RK4_KdV(Ni, c_lst(c), x0, x_limits, end_time, saving_hist, conservation, dealiasing);
        
        finalL2norm(n,c) = sqrt(abs(trapz((u_hist(:,end) - u_ana(:,end)).^2, x)));
        err2(n,c) = sum(errors.norm2)/length(errors.norm2);
        errinf(n,c) = errors.normInf(1);
        errl2(n,c) = sum(errors.L2norm)/length(errors.L2norm);
    end
end

%errl2 = errl2(:,:)./errl2(1,:);
figure('Name','Convergence plot with different velocities')
loglog(Ni_lst,errl2','-o')
legend(leg, 'Location', 'Best')
grid on
xlabel('$N$')
ylabel('$||u-\mathcal{I}_Nu||_{L2}$')


%% Question (d) - error for different domain sizes

c = 5;
x0 = 0;
d_lst = [1,2,3,4]*pi;

x = linspace(0,d_lst(end),1000);
soli = 0.5*c*sech(0.5.*sqrt(c).*(x'-x0)).^2;
figure('Name','Solition solution for domain size study')
semilogy(x, soli, 'k')
hold on
for i=1:length(d_lst)
    plot([d_lst(i), d_lst(i)], [min(soli), max(soli)], '--')
end
leg = {'Soliton', 'd = $\pi$', 'd = $2\pi$', 'd = $3\pi$', 'd = $4\pi$'};
legend(leg, 'Location', 'Best')
xlabel('Position $x$')
grid on
ylabel('Soliton $u(x,t=0)$')

saving_hist = 1;
conservation = 1;
dealiasing = 0;

end_time = 0.2;

Ni_lst = 2.^(3:8)+1;
finalL2norm = zeros(length(d_lst),length(Ni_lst));
for n=1:length(Ni_lst)
    Ni = Ni_lst(n);
    txt = ['=> Computing for N = ',num2str(Ni)];disp(txt)
    parfor d=1:length(d_lst)
        x_limits = [-d_lst(d), d_lst(d)];
        [~, ~, errors, ~, ~, ~] = RK4_KdV(Ni, c, x0, x_limits, end_time, saving_hist, conservation, dealiasing);
        finalL2norm(d, n) = sum(errors.L2norm)/length(errors.L2norm);
    end
end


cols = {[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880]};
figure('Name', 'Convergence plot with different range')
for i=1:length(d_lst)
    loglog(Ni_lst, finalL2norm(i,:), '-o', 'Color', cols{i})
    hold on
end
xlabel('N')
grid on
ylabel('$||u-\mathcal{I}_Nu||_{L2}$')
legend({'d = $\pi$', 'd = $2\pi$', 'd = $3\pi$', 'd = $4\pi$'},'Location', 'Southwest')

%% Question (e)

saving_hist = 1;
conservation = 1;
dealiasing_mode = [0,1];

Ni = 31;
c = 5;
x0 = 1;
x_limits = [-2*pi, 2*pi];
end_time = 10;
for j=1:2
    dealiasing = dealiasing_mode(j);
    [u_hist, u_ana, errors, x, time, ~] = RK4_KdV(Ni, c, x0, x_limits, end_time, saving_hist, conservation, dealiasing);

    if dealiasing
        cn_deal = fftshift(fft(u_hist));
        cn_deal = cn_deal(round((length(u_hist(:,1))-Ni)/2,0):round((length(u_hist(:,1))-Ni)/2,0)+Ni,:);
        cn_deal = cn_deal((Ni+1)/2+1:end,:);
        
        errors_deal = errors;
        
        figure
        for i=1:length(k)
            semilogy(time, (abs(cn_deal(i,:))),'DisplayName', num2str(k(i)))
            hold on
        end

    else
        errors_alias = errors;
        cn = fft(u_hist)/(Ni+1);
        cn = fftshift(cn); %Re-order accordingly
        k = 0:(Ni+1)/2-1;
        cn = cn((Ni+1)/2+1:end,:);
        
        figure
        for i=1:length(k)
            semilogy(time, (abs(cn(i,:))),'DisplayName', num2str(k(i)))
            hold on
        end

    end


    xlabel('Time')
    ylabel('$|c_n|$')
    %legend
end

figure('Name','Effect of aliasing on final solution')
semilogy(time, errors_alias.L2norm,  'DisplayName', 'Aliased')
hold on
plot(time, errors_deal.L2norm, 'DisplayName', 'Dealiased')
ylabel('$||u-\mathcal{I}_Nu||_{L2}$')
xlabel('t')
grid on
legend 

figure('Name','Std over mean')
semilogy(k,std(abs(cn'))./abs(mean(cn')), '-d')
hold on
semilogy(k,std(abs(cn_deal'))./abs(mean(cn_deal')), '-d')
legend('Baseline',"Orszag's rule")
xlabel('$n$')
ylabel('$\sigma_{\bar{c_n}}/\bar{c_n}$')
grid on


%% Question  (f)

saving_hist = 1;
conservation = 0;
dealiasing = 1;

Ni = 51;

c_arr = [0.5, 0.25];
x0_arr = [-40, -15];
x_limits = [-50, 30];
end_time = 120;

% Calling function to simulate the collision of solitons
[u_hist, errors, x, time, quant] = RK4_KdV_collision(Ni, c_arr, x0_arr, x_limits, end_time, saving_hist, conservation, dealiasing);

step = 6000;
[X,T] = meshgrid(x,time(1:step:end));
figure('Name', 'Spectral method')
surf(X,T,u_hist(:,1:step:end)')
xlabel('x')
ylabel('t')
zlabel('$\overline{u}(x,t)$')

% figure
% hold on
% for i = 1:1000:length(time)
%     plot(x,u_hist(:,i), '-r')
%     hold on
%     pause(0.001)
%     clf
% end

%% Question (g)

saving_hist = 1;
conservation = 0;

Ni_lst = 2.^(3:8) + 1;
time_arr = zeros(length(Ni_lst),2);
time_arr_norm = zeros(length(Ni_lst),2);
for j=1:2
    if j==1
        dealiasing = 0;
    else
        dealiasing = 1;
    end
    for n=1:length(Ni_lst)
        Ni = Ni_lst(n);
        txt = ['=> Computing for N = ',num2str(Ni)];
        disp(txt)
        c = 5;
        x0 = 0;
        x_limits = [-2*pi, 2*pi];
        end_time = 0.1;

        tic;
        [u_hist, u_ana, errors, x, time, quant] = RK4_KdV(Ni, c, x0, x_limits, end_time, saving_hist, conservation, dealiasing);
        time_arr(n,j) = toc; % CPU Time
        time_arr_norm(n,j) = toc/(time(2)-time(1));

    end
end
figure()
loglog(Ni_lst,time_arr(:,1), '-g', 'DisplayName', 'CPU time [ms] (aliased)')
hold on
loglog(Ni_lst, time_arr_norm(:,1), '-b', 'DisplayName', 'Normalised CPU time (aliased)')
loglog(Ni_lst,time_arr(:,2), '--g', 'DisplayName', 'CPU time [ms] (dealiased)')
loglog(Ni_lst, time_arr_norm(:,2), '--b', 'DisplayName', 'Normalised CPU time (dealiased)')
plot(Ni_lst, Ni_lst.*log(Ni_lst), '--k', 'DisplayName', '$N\log(N)$')
plot(Ni_lst, Ni_lst.^2, '-.k', 'DisplayName', '$N^2$')
grid on
xlabel('N');
ylabel('Performance');
legend('Location', 'Best')



