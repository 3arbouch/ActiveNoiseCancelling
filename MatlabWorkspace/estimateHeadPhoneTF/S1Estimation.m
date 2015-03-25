close all 
x = noiseData ; 
d = recordings(:,1); 

figure ;
plot(d) ;
title ('Recording data') ;

[STFT1,F1,T1] = spectrogram(d,1024,512,1024,song.SampleRate);
figure ; 
imagesc(T1,F1,log10(abs(STFT1)))
title('MAgnitude STFT of recordings ') ;

% x = x(44100*2:end) ; 
d = d(44100+3200:end) ; 



signalLength = length(x) ; 
filterSize= 1024; 
blockSize = filterSize ; 
numberOfIterations = 44100*3 ;  
%% NLMS
[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS, w_NLMS]=NLMS(x,d, filterSize, numberOfIterations)  ;
%% BLMS
blockSize = filterSize; 
[errorBLMS, MSerrorBLMS, timeOfConvergenceNLMSBLMS, w_BLMS]=BLMS(x,d, filterSize,blockSize, numberOfIterations)  ;
%% FDAF
numberOfIterations = numberOfIterations/filterSize -1; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterations)  ;

%% PLots



figure ; 
semilogy(MSerrorBLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS, 'r') ; 
hold on 
semilogy(MSerrorFDAF, 'k') ; 

title ('Learning curces for NLMS and BLMS') ; 
legend ('Learning curve of BLMS','Learning curve of NLMS', 'Learning curve of FDAF') ; 
ylabel('MSE') ; 
xlabel('Number of samples processed') ;



figure  ; 
bar([timeOfConvergenceNLMS, timeOfConvergenceNLMSBLMS, timeOfConvergenceNLMSFDAF]) ; 
title ('Time of computations for different Algorthims')


figure ; 
subplot(1,3,1) ; 
stem(w_NLMS) ; 
title('Transfer function estimated with NLMS') ; 
subplot(1,3,2); 
stem(w_BLMS) ; 
title('Transfer function estimated with BLMS') ; 
subplot(1,3,3) ; 
stem(w_FDAF) ; 
title('Transfer function estimated with FDAF') ; 