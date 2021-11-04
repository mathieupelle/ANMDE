function u_t = KdV(u, s,dealiasing)
    u = u';
    u_fft = fft(u);
    N = length(u);
    k = [0:N/2-1, - N/2:-1];
    u_x = real(ifft(1i*k.*u_fft))/s;
    u_xxx = real(ifft(-1i*k.^3.*u_fft))/s^3;
    u_t = -6*u.*u_x - u_xxx;
    
    
    if dealiasing
%         %Calculate original N input by the user
%         n = floor(2*N/3);
%         if mod(n,2) == 0
%             n = n-1;
%         end
% 
%         u_t_fft = fft(u_t);
%         u_t_fft = fftshift(u_t_fft);
%         u_t_fft_deal = [zeros(1,round((N-n)/2-1,0)), u_t_fft(round((N-n)/2,0):round((N-n)/2,0)+n), zeros(1,round((N-n)/2,0)-1)];
%         u_t = ifft(ifftshift(u_t_fft_deal));

        u_dea = ApplyDealiasing(u);
        u_x_dea = ApplyDealiasing(u_x);
        
        
        
        u_t = -6*u_dea.*u_x_dea - u_xxx;

       % nonlin = ifft(Dealias_Allan(fft(u),fft(u_x)));
       % u_t = -6*nonlin - u_xxx;
    end
    u_t = u_t'; 
end


