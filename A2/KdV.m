function u_t = KdV(u, s)
    u = u';
    u_fft = fft(u);
    N = length(u);
    k = [0:N/2-1, - N/2:-1];
    u_x = real(ifft(1i*k.*u_fft))/s;
    u_xxx = real(ifft(-1i*k.^3.*u_fft))/s^3;
    u_t = -6*u.*u_x - u_xxx;
    u_t = u_t';
end