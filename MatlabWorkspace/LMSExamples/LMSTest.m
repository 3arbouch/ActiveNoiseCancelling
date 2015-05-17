close all
profile on 



convergenceThreshold = 10^-4;


%%  Generate signals of interest 
signalLength = 5000 ; 
originalFilterSize = 150 ; 
[x,d,h] = generateSignals(signalLength, originalFilterSize,0 ) ; 

signalLength = 5000 ; 
originalFilterSize1 = 100 ; 
[x1,d1,h1] = generateSignals(signalLength, originalFilterSize, 1) ; 
 
figure ;
stem(h) ; 
title ('Filter Impulse reponse') ;



%% Parameters of the algorithms
addedSamples = 0 ; 
filterSize = originalFilterSize + addedSamples ; 
filterSize1 = originalFilterSize1 + addedSamples ; 
numberOfIterations = signalLength - filterSize ; 

%% NLMS 
[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS,timeOfComputationNLMS,  w_NLMS]=NLMS(x,d, filterSize, numberOfIterations, convergenceThreshold)  ;


% [errorNLMS1, MSerrorNLMS1, timeOfConvergenceNLMS1, w_NLMS1]=NLMS(x1,d1, filterSize1, numberOfIterations)  ;




% hold on 
% semilogy(MSerrorNLMS1, 'r') ; 
% legend('Filter Lenght= 100', 'Filter Lenght= 200')
% title ('Learning curve  of NLMS for different filter Sizes') ;
% ylabel('MSE') ; 
% xlabel('Number of samples processed') ;
% 
% figure ; 
% plot(errorNLMS, 'b') ; 
% title('NLMS Error for filter Size 100') ; 
% ylabel('Error') ; 
% xlabel('Number of iterations') ;





%% Block LMS
blockSize = filterSize; 
[errorBLMS, MSerrorBLMS, timeOfConvergenceNLMSBLMS,timeOfComputationBLMS,  w_BLMS]=BLMS(x,d, filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;



%% FDAF
blockSize = filterSize/2 ; 
numberOfIterations = signalLength/filterSize -1 ; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;
[errorFDAFCC, MSerrorFDAFCC, timeOfConvergenceNLMSFDAFCC,timeOfComputationFDAFCC, w_FDAFCC]=FDAFCC(x,d, filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;


figure ;
plot(errorFDAF);

figure  ; 
bar([timeOfConvergenceNLMS timeOfComputationNLMS ; timeOfConvergenceNLMSBLMS timeOfComputationBLMS ; timeOfConvergenceNLMSFDAF timeOfComputationFDAF ; timeOfConvergenceNLMSFDAFCC timeOfComputationFDAFCC]) ; 
title ('Time of convergence and computation for different Algorthims')



figure ; 
semilogy(MSerrorBLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS, 'r') ; 
hold on 
semilogy(MSerrorFDAF, 'k') ; 
hold on 
semilogy(MSerrorFDAFCC, 'g') ; 

title ('Learning curces for NLMS and BLMS with filter Size= 100') ; 
legend ('Learning curve of BLMS','Learning curve of NLMS', 'Learning curve of FDAF') ; 
ylabel('MSE') ; 
xlabel('Number of samples processed') ;



% PLot the estimated transferFunction with different methods

figure ; 
subplot(3,2,1) ; 
stem(h) ; 
title('Original Transfer function') ; 
subplot(3,2,2) ; 
stem(w_NLMS) ; 
title('Transfer function estimated with NLMS') ; 

subplot(3,2,3); 
stem(w_BLMS) ; 
title('Transfer function estimated with BLMS') ; 
subplot(3,2,4) ; 
stem(w_FDAF) ; 
title('Transfer function estimated with FDAF') ; 
subplot(3,2,5) ; 
stem(w_FDAFCC) ; 
title('Transfer function estimated with FDAFCC') ; 


figure ;


profile viewer

%%
figure ; 
plot(filter(w_FDAF,1,x));
hold on ;
plot(d);