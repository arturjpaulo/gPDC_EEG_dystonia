
function plot_periodogram_fft (raw)

N = length(raw.data);
xdft = fft(raw.data);
xdft = xdft(1:N/2+1);
psdx = (1/(raw.Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:raw.Fs/length(raw.data):raw.Fs/2;

plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

end