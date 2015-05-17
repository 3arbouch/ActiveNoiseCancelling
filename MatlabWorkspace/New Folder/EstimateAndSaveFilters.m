%% This script is used to estimate and save the filters from the saved experiments
clear all 
close all
data = load('experiment30.mat');


experimentName = 'experiment30Filters(4096)';
% Size of the estimated filters
filterSize = 1024*4 ;
blockSize = filterSize ;
convergenceThreshold = 10^-6 ;


%% Estimation of s1(z)
referenceSignal = data.noise' ;
filtredSignal = data.referenceMic';

numberOfIterations = length(referenceSignal)/filterSize -1 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  S1 ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold );



%% PLots

figure ;
stem(S1);
title ('Estimated S1(z)')
figure ;
semilogy(MSerror);
title ('Mean square error')

echoCancellerSignal = filter(S1,1,referenceSignal);
echoCancelledSignal = filtredSignal(1:length(referenceSignal)) - echoCancellerSignal ;
figure ;
plot(filtredSignal);
hold on ;
plot(echoCancellerSignal);



%% Estimation of S2(z)
referenceSignal = data.noise' ;
filtredSignal = data.errorMic';

numberOfIterations = length(referenceSignal)/filterSize -1 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  S2 ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold );


 %% Remove the delay from S2(z)
 [~, index]= max(S2);
 S2 = S2(index:end);
%% Plots
figure ;
stem(S2);
title ('Estimated S2(z)')
figure ;
semilogy(MSerror);
title ('Mean square error')


estimatedFiltredSignal = filter(S2, 1, referenceSignal);
figure; 
plot(filtredSignal);
hold on 
plot(estimatedFiltredSignal);

%% Estimation of P(z)
referenceSignal = data.referenceMicPrime';
filtredSignal = data.errorMicPrime';
referenceSignal_prime= filter(S2,1,referenceSignal);

numberOfIterations = length(referenceSignal)/filterSize -1 ;
% numberOfIterations = 50 ;
[ error, MSerror, timeOfConvergence,timeOfComputation,  P ] = FDAFOSM( referenceSignal_prime, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold );


figure ;
stem(P);
title ('Estimated P(z)')
figure ;
semilogy(MSerror);
title ('Mean square error')


estimatedFiltredSignal = filter(P,1,referenceSignal_prime);
figure ;
plot(filtredSignal);
hold on 
plot(estimatedFiltredSignal);


%% Save the filters 

save(experimentName,'S1','S2','P');
