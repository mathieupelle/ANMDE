%% Advanced Numerical Methods for Differential Equations - Assignment 2
clear variables; close all; clc; beep off;
%% Question (a) - Legendre Tau Method (LTM)

%Input: 
N = [2 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9]; %how many basis functions?
eps = [0.1 0.01 0.001]; %Function parameter
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

%% Question (a) - Legendre Collocation Method


