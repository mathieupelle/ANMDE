function [u_hist, u_ana, errors, x, time, quant] = RK4_KdV(Ni, c, x0, saving_hist, conservation)
    
    % Spatial discretisation
    x_end = 2*pi; %left xlim
    x_start = -2*pi; %right xlim
    scaling = (x_end-x_start)/(2*pi); %scaling factor
    dx = (x_end-x_start)/(Ni+1); %spatial element width
    x = x_start:dx:x_end-dx; %x vector
    
    % Temporal discretisation and IC
    u_init = 0.5*c*sech(0.5*sqrt(c).*(x-c*0-x0)).^2; %initial guess
    dt = 2.82/(3*(Ni+1)*max(abs(u_init))+(Ni+1)^3/8); %time step from RK4 stability
    time = 0:dt:1; %time 
    
    
    % Allocating storage
    [u_hist, u_ana] = deal(zeros(length(x), length(time))); %spectral and analytic solutions  
    u = u_init';
    u_hist(:,1) = u_init;
    u_ana(:,1) = u_init;
    
    [L2norm, norm2, normInf, M, M_ana, V, V_ana, E, E_ana] = deal(zeros(length(time), 1));
    
    L2norm(1) = sqrt(trapz((u_hist(:,1) - u_ana(:,1)).^2, x));
    
    if conservation == 1
        % Errors
        norm2(1,1) = norm(u - u_ana(:,1) ,2); %2 Norm
        normInf(1,1) = norm(u - u_ana(:,1) , Inf); %Infinity Norm
        
        %Conservation quantities
        M(1,1) = trapz(u, x); %Mass (spectral)
        M_ana(1,1) = trapz(u_ana(:,1), x); %Mass (analytical)
        V(1,1) = trapz(u.^2, x); %Momentum (spectral)
        V_ana(1,1) = trapz(u_ana(:,1).^2, x); %Momentum (analytical)
        D = get_FourierDiffMatrix(Ni+1,false)/scaling; %Differentiation matrix
        ux = D*u; %Derivative (spectral)
        ux_ana = D*u_ana(:,1); %Derivative (analytical)
        E(1,1) = trapz(0.5*ux.^2 - u.^3, x); %Momentum (spectral)
        E_ana(1,1) = trapz(0.5*ux_ana.^2 - u.^3, x); %Momentum (analytical)

    end

%     D = get_FourierDiffMatrix(Ni+1,false)/2;
%     D3 = D*D*D;
    for t=1:length(time)-1

        % Diff. matrix Runge-Kutta
%         K1 = -6*u'*D*u - D3*u;
%         K2 = -6*(u'+dt*K1/2)*D*(u+dt*K1/2) - D3*(u+dt*K1/2);
%         K3 = -6*(u'+dt*K2/2)*D*(u+dt*K2/2) - D3*(u+dt*K2/2);
%         K4 = -6*(u'+dt*K3/2)*D*(u+dt*K3/2) - D3*(u+dt*K3/2);

        %FFT Runge-Kutta
        K1 = KdV(u, scaling);
        K2 = KdV(u+dt*K1/2, scaling);
        K3 = KdV(u+dt*K2/2, scaling);
        K4 = KdV(u+dt*K3, scaling);
        u = u + dt/6*(K1+2*K2+2*K3+K4);
        
        if saving_hist == 1 %saving each time step
            u_hist(:,t+1) = u;
            x_ana = x+c*time(t);
            
            d = x_end-x_start;

            %u_ana(:,t+1) = 0.5*c*sech(0.5*sqrt(c)*(x-c*time(t+1)-x0)).^2; %analytic
            u_ana(:,t+1) = 0.5*c*sech(0.5*sqrt(c)*(x-c*time(t+1)-x0)).^2 ...
            + 0.5*c*sech(0.5*sqrt(c)*(x-c*time(t+1)-x0+d)).^2 ...
            + 0.5*c*sech(0.5*sqrt(c)*(x-c*time(t+1)-x0+2*d)).^2 ... 
            + 0.5*c*sech(0.5*sqrt(c)*(x-c*time(t+1)-x0+3*d)).^2; %analytic

%             for k=0:100
%                 if (x_end-dx)*(k+1) > x_ana(end)
%                     break
%                 end
%             end
%             %Calculate necessary distances for the reassignemnt
%             l = x_ana(end)-(x_end-dx)*k;
%             d = ((x_end-dx) - x_start)/2;
%             
%             
%             u_ana(:,t+1) = [u_ana(x_ana>(x_ana(end)-l),t+1); u_ana(x_ana<=(x_ana(end)-l),t+1)];       
%             
            
            L2norm(t+1) = sqrt(abs(trapz((u - u_ana(:,t+1)).^2, x)));  %L2 norm
            
            if conservation == 1
                norm2(t+1) = norm(u - u_ana(:,t+1) ,2);
                normInf(t+1) = norm(u - u_ana(:,t+1) , Inf);
                M(t+1) = trapz(u, x);
                M_ana(t+1) = trapz(u_ana(:,t+1), x);
                V(t+1) = trapz(u.^2, x);
                V_ana(t+1) = trapz(u_ana(:,t+1).^2, x);
                ux = D*u;
                ux_ana = D*u_ana(:,t+1);
                E(t+1) = trapz(0.5*ux.^2 - u.^3, x);
                E_ana(t+1) = trapz(0.5*ux_ana.^2 - u.^3, x);
            end
        else 
            if t == length(time)-1 %saving last time step only
                u_hist(:,end) = u; %spectral
                u_ana(:,end) = 0.5*c*sech(0.5*sqrt(c)*(x-c*time(end)-x0)).^2; %analytic
            end
        end

    end
    
    errors.L2norm = L2norm;
    errors.norm2 = norm2;
    errors.normInf = normInf;
    quant.M = M;
    quant.V = V;
    quant.E = E;
    quant.M_ana = M_ana;
    quant.V_ana = V_ana;
    quant.E_ana = E_ana;
    
end  