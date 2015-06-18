%% This script is used to estimate and save the filters from the saved experiments
clear all 
close all
data = load('experiment100(8Khz).mat');


experimentName = 'experiment100Filters(1024, 8khz)';
% Size of the estimated filters
filterSize = 1024 ;
blockSize = filterSize ;
convergenceThreshold = 10^-11 ;

fs = 8000 ;
%% Estimation of s1(z)
referenceSignal = data.noise' ;
filtredSignal = data.referenceMic';

numberOfIterations = length(referenceSignal)/filterSize -1 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  S1 ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,  numberOfIterations, convergenceThreshold );



%% PLots



tfplot(S1,fs,'S1','Estimated echo canceller filter S1(z)')
figure ;
semilogy(MSerror);
title ('Mean square error of S1(z)')

echoCancellerSignal = filter(S1,1,referenceSignal);
echoCancelledSignal = filtredSignal(1:length(referenceSignal)) - echoCancellerSignal ;
figure ;
plot(filtredSignal);
hold on ;
plot(echoCancellerSignal);



%% Estimation of S2(z)
referenceSignal = data.noise' ;
filtredSignal = data.errorMic';

numberOfIterations = length(referenceSignal)/filterSize -2 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  S2 ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,  numberOfIterations, convergenceThreshold );


 

%% Plots

tfplot(S2,fs,'S2','Estimated inside headphone filter S2(z)');

figure ;
semilogy(MSerror);
title ('Mean square error of S2(z)')


estimatedFiltredSignal = filter(S2, 1, referenceSignal);
figure; 
plot(filtredSignal);
hold on 
plot(estimatedFiltredSignal);

%% Remove the delay from S2(z)
 [~, index]= max(S2);
 S2Prime = S2 ; 
  S2 = S2(index:end);
%% Estimation of P(z)
referenceSignal = data.referenceMicPrime';
filtredSignal = data.errorMicPrime';
referenceSignal_prime= filter(S2,1,referenceSignal);

numberOfIterations = length(referenceSignal)/filterSize -5 ;
% numberOfIterations = 50 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  P ] = FDAFOSM( referenceSignal_prime, filtredSignal, filterSize,  numberOfIterations, convergenceThreshold );


tfplot(P,fs,'P','Estimated headphone filter P(z)')

figure ;
semilogy(MSerror);
title ('Mean square error of P(z)')


estimatedFiltredSignal = filter(P,1,referenceSignal_prime);
figure ;
plot(filtredSignal);
hold on 
plot(estimatedFiltredSignal);


%% Save the filters 

save(experimentName,'S1','S2','P', 'S2Prime');
