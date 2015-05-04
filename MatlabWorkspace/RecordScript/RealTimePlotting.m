Fs = 8000;                    %# sampling frequency in Hz
T = 0.01;                        %# length of one interval signal in sec
t = 0:1/Fs:T-1/Fs;            %# time vector
nfft = 2^nextpow2(Fs);        %# n-point DFT
numUniq = ceil((nfft+1)/2);   %# half point
f = (0:numUniq-1)'*Fs/nfft;   %'# frequency vector (one sided)

%# prepare plots
figure
hAx(1) = subplot(211);
hLine(1) = line('XData',t, 'YData',nan(size(t)), 'Color','b', 'Parent',hAx(1));
xlabel('Time (s)'), ylabel('Amplitude')
hAx(2) = subplot(212);
hLine(2) = line('XData',f, 'YData',nan(size(f)), 'Color','b', 'Parent',hAx(2));
xlabel('Frequency (Hz)'), ylabel('Magnitude (dB)')
set(hAx, 'Box','on', 'XGrid','on', 'YGrid','on')
%#specgram(sig, nfft, Fs);

%# prepare audio recording
recObj = audiorecorder(Fs,8,1);

%# Record for 10 intervals of 1sec each
disp('Start speaking...')
for i=1:1000
    recordblocking(recObj, T);

    %# get data and compute FFT
    sig = getaudiodata(recObj);
    fftMag = 20*log10( abs(fft(sig,nfft)) );

    %# update plots
    set(hLine(1), 'YData',sig)
    set(hLine(2), 'YData',fftMag(1:numUniq))
    title(hAx(1), num2str(i,'Interval = %d'))
    drawnow                   %# force MATLAB to flush any queued displays
end
disp('Done.')
