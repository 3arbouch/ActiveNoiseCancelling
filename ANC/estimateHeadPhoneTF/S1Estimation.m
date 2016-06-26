close all 
profile on

load('Workspace.mat');
x = noiseData ; 
d = recordings(:,1) ;

convergenceThreshold = 10^-4;

figure ;
plot(d) ;
title ('Recording data') ;

[STFT1,F1,T1] = spectrogram(d,1024,512,1024,song.SampleRate);
figure ; 
imagesc(T1,F1,log10(abs(STFT1)))
title('MAgnitude STFT of recordings ') ;
% Compute the corrolation 
datacorr = xcorr(d,noiseData);
[value, index]= max(datacorr) ;
offset = 20 ;
startPoint = index -length(noiseData)- offset;

d = d(startPoint:end) ; 



  
%% NLMS

signalLength = length(x) ; 
filterSize= 200; 
blockSize = filterSize ; 
numberOfIterations = 44100*5 ;

[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS,timeOfComputationNLMS,  w_NLMS]=NLMS(x,d, filterSize, numberOfIterations, convergenceThreshold)  ;

%% BLMS
blockSize = filterSize; 
[errorBLMS, MSerrorBLMS, timeOfConvergenceNLMSBLMS,timeOfComputationBLMS,  w_BLMS]=BLMS(x,d, filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;
% figure ; 
% semilogy(MSerrorBLMS, 'b') ;
% title ('NLMS learning curve'); 
% 
% 
% figure ; 
% stem(abs(w_BLMS)); 
% title ('Estimated fiter with NLMS');

%% FDAF
numberOfIterationsF = numberOfIterations/filterSize -1; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterationsF, convergenceThreshold)  ;
[errorFDAFCC, MSerrorFDAFCC, timeOfConvergenceNLMSFDAFCC,timeOfComputationFDAFCC, w_FDAFCC]=FDAFCC(x,d, filterSize,blockSize, numberOfIterationsF, convergenceThreshold)  ;



%% PLots

figure ;
plot(errorFDAF);

figure ; 
semilogy(MSerrorBLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS, 'r') ; 
hold on 
semilogy(MSerrorFDAF, 'k') ; 
hold on 
semilogy(MSerrorFDAFCC, 'g') ; 
hold on 
y = convergenceThreshold*ones(numberOfIterations,1);
semilogy(y, 'm');


title ('Learning curces for NLMS and BLMS') ; 
legend ('Learning curve of BLMS','Learning curve of NLMS', 'Learning curve of FDAF', 'Learning curve of FDAFCC') ; 
ylabel('MSE') ; 
xlabel('Number of samples processed') ;



figure  ; 
bar([ timeOfComputationNLMS ;  timeOfComputationBLMS ;  timeOfComputationFDAF ;  timeOfComputationFDAFCC], 0.5, 'r') ; 
title ('Time of computation  for different Algorthims')

figure  ; 
bar([timeOfConvergenceNLMS  ; timeOfConvergenceNLMSBLMS  ; timeOfConvergenceNLMSFDAF  ; timeOfConvergenceNLMSFDAFCC ], 0.5, 'y') ; 
title ('Time of convergence  for different Algorthims')

figure ; 
subplot(2,2,1) ; 
stem(w_NLMS) ; 
title('Transfer function estimated with NLMS') ; 
subplot(2,2,2); 
stem(w_BLMS) ; 
title('Transfer function estimated with BLMS') ; 
subplot(2,2,3) ; 
stem(w_FDAF) ; 
title('Transfer function estimated with FDAF') ; 

subplot(2,2,4) ; 
stem(w_FDAFCC) ; 
title('Transfer function estimated with FDAFCC') ; 

profile off 
profile viewer
