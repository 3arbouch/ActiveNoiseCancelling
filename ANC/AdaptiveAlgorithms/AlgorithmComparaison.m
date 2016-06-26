% This script is used to perform the a comparaison in terms of the
% convergance rate, the convergence time and the computation time 
% between four different
% adaptive algorithms: NLMS, BLMS , FDAFOSM and FDAFCC
clear all 
close all 
profile on 





%%  Generate signals of interest 
signalLength = 25000 ; 
originalFilterSize = 128 ; 
[x,d,h] =  generateSignals( signalLength, originalFilterSize, 'wn', 'ma', 0 );

signalLength = 25000 ; 
originalFilterSize1 = 512 ; 
[x1,d1,h1] = generateSignals( signalLength, originalFilterSize, 'wn', 'ma', 1 );

 
figure ;
stem(h) ; 
title ('Filter Impulse reponse') ;



%% Parameters of the algorithms
addedSamples = 0 ; 
filterSize = originalFilterSize + addedSamples ; 
filterSize1 = originalFilterSize1 + addedSamples ; 
numberOfIterations = signalLength - filterSize ; 
convergenceThreshold = 10^-5;

%% NLMS 
[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS,timeOfComputationNLMS,  w_NLMS]=NLMS(x,d, filterSize, numberOfIterations, convergenceThreshold)  ;

% This portion of code compares the performance of the algorithm for two
% different filter sizes
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
numberOfIterations = signalLength/filterSize -1 ; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(x,d, filterSize, numberOfIterations, convergenceThreshold)  ;
[errorFDAFCC, MSerrorFDAFCC, timeOfConvergenceNLMSFDAFCC,timeOfComputationFDAFCC, w_FDAFCC]=FDAFCC(x,d, filterSize, numberOfIterations, convergenceThreshold)  ;


 


%% Plot  bar plots 
% Plot a hitogram showing the corgence and the computation time of each
% algorithm 


figure; 
timeOfconvergences = [timeOfConvergenceNLMS  ; timeOfConvergenceNLMSBLMS  ; timeOfConvergenceNLMSFDAF   ] ; 
b = bar(timeOfconvergences,'FaceColor','r','LineWidth',2) ; 
xlabel('Algorithms');
ylabel('Seconds');
title ('Convergence Time')

Labels = {'NLMS', 'BLMS', 'FDAFOSM'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels, 'fontsize', 30);
 yb = cat(1, b.YData);
 xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
 hold on;
 text(xb(:),yb(:), cellstr(num2str(timeOfconvergences(:))), 'fontsize',30);



figure; 
timeOfcomputions = [timeOfComputationNLMS  ; timeOfComputationBLMS  ; timeOfComputationFDAF   ] ; 
b = bar(timeOfcomputions,'FaceColor','y','LineWidth',2) ; 
xlabel('Algorithms');
ylabel('Seconds');
title('Computation time')

Labels = {'NLMS', 'BLMS', 'FDAFOSM'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels, 'fontsize', 30);
 yb = cat(1, b.YData);
 xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
 hold on;
 text(xb(:),yb(:), cellstr(num2str(timeOfcomputions(:))), 'fontsize',30);


 
%% Plot learning curves 

figure ; 
semilogy(MSerrorBLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS, 'r') ; 
hold on 
semilogy(MSerrorFDAF, 'k') ; 


 title ('Learning curves for NLMS, BLMS and FDAF with filter Size= 150') ; 
legend ('Learning curve of BLMS','Learning curve of NLMS', 'Learning curve of FDAF') ; 
ylabel('MSE') ; 
xlabel('Number of samples processed') ;


%%  PLot the estimated transferFunction with different methods

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




profile viewer

