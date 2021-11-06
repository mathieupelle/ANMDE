function what=Dealias_Allan(uhat,vhat)
N=length(uhat);
M=length(vhat);
uhatpad = [uhat(1:N/2) zeros(1, M-N) uhat(N/2+1:end)];
vhatpad = [vhat(1:N/2) zeros(1, M-N) vhat(N/2+1:end)];
upad = ifft(uhatpad);
vpad = ifft(vhatpad);
wpad = upad.*vpad;
wpad_hat = fft(wpad);
what = 3/2*[wpad_hat(1:N/2) wpad_hat(M-N/2+1:M)];