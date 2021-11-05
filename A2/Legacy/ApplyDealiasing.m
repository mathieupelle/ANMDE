function x_dealiased = ApplyDealiasing(x)

N = length(x);

n = floor(2*N/3);
if mod(n,2) == 0
    n = n-1;
end

x_fft = fft(x);
x_fft = fftshift(x_fft);
x_fft_deal = [zeros(1,round((N-n)/2-1,0)), x_fft(round((N-n)/2,0):round((N-n)/2,0)+n), zeros(1,round((N-n)/2,0)-1)];
x_dealiased = ifft(ifftshift(x_fft_deal));

end