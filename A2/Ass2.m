%% Advanced Numerical Methods for Differential Equations - Assignment 2
clear variables; close all; clc; beep off;

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','normal')
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
conservation = 0;

Ni_lst = 2.^(3:9) + 1;
Ni_lst = 31;
finalL2norm = zeros(length(Ni_lst),1);
for n=1:length(Ni_lst)
    Ni = Ni_lst(n);
    txt = ['=> Computing for N = ',num2str(Ni)];
    disp(txt)
    c = 5;
    x0 = 0;

    [u_hist, u_ana, errors, x, time, quant] = RK4_KdV(Ni, c, x0, saving_hist, conservation);
    
    finalL2norm(n,1 ) = sqrt(abs(trapz((u_hist(:,end) - u_ana(:,end)).^2, x)));

    if n==length(Ni_lst)
        if saving_hist ==1

%         [X,T] = meshgrid(x,time);
%         figure('Name', 'Spectral method')
%         surf(X,T,u_hist')
%         xlabel('x')
%         ylabel('t')
%         zlabel('$\overline{u}(x,t)$')
%         
%         figure('Name', 'Soliton')
%         surf(X,T,u_ana')
%         xlabel('x')
%         ylabel('t')
%         zlabel('u(x,t)')

%         figure('Name', 'Absolute error')
%         surf(X,T,abs(u_ana'-u_hist'))
%         xlabel('x')
%         ylabel('t')
%         zlabel('$|\overline{u}(x,t) - u(x,t)|$')

            figure('Name', 'L2 norm in time')
            plot(time, errors.L2norm)
            xlabel('Time')
            ylabel('L2')
            grid on
        end
        figure('Name', 'Final timestep')
        plot(x, u_ana(:,end), 'DisplayName', 'analytic')
        hold on
        plot(x, u_hist(:,end), '--.k', 'DisplayName', 'spectral method')
        xlabel('x')
        ylabel('u')
        legend
        grid on
    end
end

figure('Name', 'Convergence')
loglog(Ni_lst, finalL2norm)
grid on
xlabel('N')
ylabel('$||u-\mathcal{I}_Nu||_{L2}$')

% if saving_hist == 1
%     figure
%     hold on
%     for i = 1:10:length(time)
%         plot(x,u_hist(:,i), '-r')
%         hold on
%         plot(x,u_ana(:,i), '--ok')
%         pause(0.005)
%         clf
%     end
% end

%% Question (d)

saving_hist = 1;
conservation = 1;

c_lst = [0.25, 0.5, 1];
leg = {'c = 0.25', 'c = 0.5', 'c = 1.0'};
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
for i=1:length(c_lst)
    
    Ni = 55;
    x0 = 0;
    c = c_lst(i);
    
    [~, ~, errors, ~, time, quant] = RK4_KdV(Ni, c, x0, saving_hist, conservation);
    
    set(0, 'CurrentFigure', f1)
    loglog(time, errors.norm2)
    grid on
    hold on
    xlabel('Time')
    ylabel('$||u-\mathcal{I}_Nu||_2$')
    legend(leg, 'Location', 'Best')
    
    set(0, 'CurrentFigure', f2)
    loglog(time, errors.normInf)
    grid on
    hold on
    xlabel('Time')
    ylabel('$||u-\mathcal{I}_Nu||_\infty$')
    legend(leg, 'Location', 'Best')

    set(0, 'CurrentFigure', f3)
    err = abs(quant.M - quant.M_ana)./abs(quant.M_ana);
    loglog(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Mass (M)')
    legend(leg, 'Location', 'Best')
    
    set(0, 'CurrentFigure', f4)
    err = abs(quant.V - quant.V_ana)./abs(quant.V_ana);
    loglog(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Momentum (V)')
    legend(leg, 'Location', 'Best')
    
    
    set(0, 'CurrentFigure', f5)
    err = abs(quant.E - quant.E_ana)./abs(quant.E_ana);
    loglog(time, err)
    grid on
    hold on
    xlabel('Time')
    ylabel('Energy (E)')
    legend(leg, 'Location', 'Best')
    
end

%% Question (e)

saving_hist = 1;
conservation = 0;

Ni_lst = 3:4:51;
%col = {'r', 'g', 'k', 'b', 'm'};
figure
for i=1:length(Ni_lst)
    
    Ni = Ni_lst(i);
    c = 5;
    x0 = 0;

    [u_hist, u_ana, ~, x, time, ~] = RK4_KdV(Ni, c, x0, saving_hist, conservation);

    cn = fft(u_hist)/(Ni+1);
    cn = fftshift(cn); % Re-order accordingly
    k = 0:(Ni+1)/2-1;
    cn = cn((Ni+1)/2+1:end,:);

%     cn_ana = fft(u_ana)/(Ni+1);
%     cn_ana = fftshift(cn_ana); % Re-order accordingly
%     cn_ana = cn_ana((Ni+1)/2+1:end,:);
    k_nyq = pi/(x(2)-x(1));
    

    plot(k, sum(abs(cn),2), 'DisplayName', append('N = ', num2str(Ni)))
    hold on
    %plot([k_nyq, k_nyq], [0, 0.4])
    

end

grid on
xlabel('N')
ylabel('$c_n$')
xlim([0,15])
legend

figure
for i=1:length(k)
    plot(time, abs(cn(i,:))-mean(abs(cn(i,:))))
    hold on
end
xlabel('Time')
ylabel('$c_n$')

%% Question  (f)

saving_hist = 1;
conservation = 0;

Ni = 55;

c_arr = [0.5, 0.25];
x0_arr = [-40, -15];

% Calling function to simulate the collision of solitons
[u_hist, errors, x, time, quant] = RK4_KdV_collision(Ni, c_arr, x0_arr, saving_hist, conservation);

[X,T] = meshgrid(x,time(1:800:end));
figure('Name', 'Spectral method')
surf(X,T,u_hist(:,1:800:end)')
xlabel('x')
ylabel('t')
zlabel('$\overline{u}(x,t)$')

% figure
% hold on
% for i = 1:500:length(time)
%     plot(x,u_hist(:,i), '-r')
%     hold on
%     pause(0.1)
%     clf
% end

%% Question (g)

saving_hist = 1;
conservation = 0;

Ni_lst = 2.^(3:9) + 1;
time_arr = zeros(length(Ni_lst),1);

for n=1:length(Ni_lst)
    Ni = Ni_lst(n);
    txt = ['=> Computing for N = ',num2str(Ni)];
    disp(txt)
    c = 5;
    x0 = 0;

    tic;
    [u_hist, u_ana, errors, x, time, quant] = RK4_KdV(Ni, c, x0, saving_hist, conservation);
    time_arr(n,1) = toc; % CPU Time
    
end

%%

figure()
loglog(Ni_lst,time_arr)
grid on
xlabel('N');
ylabel('CPU Time [ms]');


