close all 
clear all 
load('Experiment1.mat');
x = recordings(:,1); 
d = recordings(:,2) ;

x = x(44100:end);
d= d(44100:end);
figure ;
plot(d) ;
title ('Recording data') ;

[STFT1,F1,T1] = spectrogram(d,1024,512,1024,44100);
figure ; 
imagesc(T1,F1,log10(abs(STFT1)))
title('MAgnitude STFT of recordings ') ;
% Compute the corrolation 
% datacorr = xcorr(d,noiseData);
% [value, index]= max(datacorr) ;
% offset = 20 ;
% startPoint = index -length(noiseData)- offset;
% 
% d = d(startPoint:end) ; 



  
%% NLMS

signalLength = length(x) ; 
filterSize= 44100; 
blockSize = filterSize ; 
numberOfIterations = 44100*5 ;

[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS, w_NLMS]=NLMS(x,d, filterSize, numberOfIterations)  ;

figure ; 
semilogy(MSerrorNLMS, 'b') ;
title ('NLMS learning curve'); 


figure ; 
stem(abs(w_NLMS)); 
title ('Estimated fiter with NLMS');

%% BLMS
blockSize = filterSize; 
[errorBLMS, MSerrorBLMS, timeOfConvergenceNLMSBLMS, w_BLMS]=BLMS(x,d, filterSize,blockSize, numberOfIterations)  ;
% figure ; 
% semilogy(MSerrorBLMS, 'b') ;
% title ('NLMS learning curve'); 
% 
% 
% figure ; 
% stem(abs(w_BLMS)); 
% title ('Estimated fiter with NLMS');

%% FDAF
numberOfIterations = numberOfIterations/filterSize -1; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterations)  ;
figure ;
semilogy(MSerrorFDAF, 'k') ; 

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