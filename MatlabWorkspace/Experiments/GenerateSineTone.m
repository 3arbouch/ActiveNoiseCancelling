close all 

f0= 8 ;
f1 = 15000 ;
fs = 44100;
sineTime = 2 ;

time = 0:(1/fs):sineTime - (1/fs) ;
sine1= sin(2*pi*f0.*time);
sine2 =  sin(2*pi*f1.*time);

signal = sine1 + sine2 ;

step = 1/length(time) ;
frequencyAxis = [(-1/2):1/length(time):(1/2) - step];
signalFFT = abs(fft(signal));
signalFFT = 1/(sineTime*fs/2).*[signalFFT(length(signalFFT)/2+1:end), signalFFT(1:length(signalFFT)/2)];

subplot(2,1,1)
plot(frequencyAxis,signalFFT);
subplot(2,1,2)
plot(signal);


sound(signal,fs);
