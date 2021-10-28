%% Advanced Numerical Methods for Differential Equations - Assignment 2
clear variables; close all; clc; beep off;
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
Ni = 2^3; 

% Create coordinates array
r1 = 1; 
a = 3;
r2 = a*r1;
r = linspace(r1,r2,Ni+1);
%theta = 0:2*pi/(Ni+1):(2*pi);
N_step = 2*pi/(Ni+1);
theta = 0:N_step:2*pi-N_step;
%theta = theta(1:end-1);

%Create meshgrid
[R,TH] = meshgrid(r,theta);

%Vectorize mesghrid
R = R(:); TH = TH(:);

%Create differentation matrices
% For radial:
x_basis = JacobiGL(0, 0, Ni); % Get the nodes (Gauss Lobato), using both endpoints (not idiot)

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
%Dr(1,:) = [1,zeros(1,Ni)];
%Dr(end,:) = [zeros(1,Ni),1];

%Get differentation matrix with Fourier basis
Dth = get_FourierDiffMatrix(Ni+1,false);
Dthth = Dth*Dth;

%Create them in 2D
I = eye(Ni+1);
%DR = kron(I,Dr); 
%DTT = kron(Dthth,I);
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
%L(1:Ni+1,1) = 1; L(Ni+2:end,1) = 0;
%L(end-Ni:end,end) = 1; L(1:end-Ni-1,end) = 0;
%f(1:Ni+1:end) = 1;
%f(Ni+1:Ni+1:end) = 1;

%Solve the system
u = L\f;
u_reshaped = reshape(u(:),Ni+1,Ni+1);

%Compute the error
[R,TH] = meshgrid(r,theta);
u_th = V_inf*(R+r1^2*R.^-1).*cos(TH);
err = abs((u_reshaped-u_th)./u_th);

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






