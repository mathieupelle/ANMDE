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

            [N,tN] = size(u');
            K = N/2*3; % 3/2 zero padding
            upad = zeros(K,tN);
            uxpad = zeros(K,tN);
            indvpad = [1:N/2,K-N/2+1:K];
            upad(indvpad,:) = u;
            uxpad(indvpad,:) = u_x;
            temp = fft(real(ifft(upad.*uxpad)));
            temp = K/N *temp; % scale back to N-FFT
            Fv2 = temp(indvpad); 
            Fv2(N/2+1) = 0; % remove the padding zeros
            nonlin = Fv2';

       %nonlin = ApplyDealiasing(u).*ApplyDealiasing(u_x);
       u_t = -6*nonlin - u_xxx;
    end
    u_t = u_t'; 
end


