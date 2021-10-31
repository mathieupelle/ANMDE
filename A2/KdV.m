function u_t = KdV(u)
    u = u';
    u_fft = fft(u);
    N = length(u);
    k = [0:round(N/2-1), - round(N/2):-1];
    u_x = real(ifft(1i*k.*u_fft));
    u_xxx = real(ifft(1i*k.^3.*u_fft));
    u_t = -6*u.*u_x + u_xxx;
    u_t = u_t';
end