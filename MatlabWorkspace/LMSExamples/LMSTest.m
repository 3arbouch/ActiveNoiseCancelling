close all
%%  Generate signals of interest 
signalLength = 5000 ; 
originalFilterSize = 100 ; 
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
[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS, w_NLMS]=NLMS(x,d, filterSize, numberOfIterations)  ;
[errorNLMS1, MSerrorNLMS1, timeOfConvergenceNLMS1, w_NLMS1]=NLMS(x1,d1, filterSize1, numberOfIterations)  ;

figure ; 
semilogy(MSerrorNLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS1, 'r') ; 
legend('Filter Lenght= 100', 'Filter Lenght= 200')
title ('Learning curve  of NLMS for different filter Sizes') ;
ylabel('MSE') ; 
xlabel('Number of samples processed') ;

figure ; 
plot(errorNLMS, 'b') ; 
title('NLMS Error for filter Size 100') ; 
ylabel('Error') ; 
xlabel('Number of iterations') ;





%% Block LMS
blockSize = filterSize; 
[errorBLMS, MSerrorBLMS, timeOfConvergenceNLMSBLMS, w_BLMS]=BLMS(x,d, filterSize,blockSize, numberOfIterations)  ;



%% FDAF
blockSize = filterSize/2 ; 
numberOfIterations = signalLength/filterSize -1 ; 
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterations)  ;

figure  ; 
plot(errorFDAF) ; 
title('Error with FDAF')


figure  ; 
bar([timeOfConvergenceNLMS, timeOfConvergenceNLMSBLMS, timeOfConvergenceNLMSFDAF]) ; 
title ('Time of convergence for different Algorthims')



figure ; 
semilogy(MSerrorBLMS, 'b') ;
hold on 
semilogy(MSerrorNLMS, 'r') ; 
hold on 
semilogy(MSerrorFDAF, 'k') ; 

title ('Learning curces for NLMS and BLMS with filter Size= 100') ; 
legend ('Learning curve of BLMS','Learning curve of NLMS', 'Learning curve of FDAF') ; 
ylabel('MSE') ; 
xlabel('Number of samples processed') ;



% PLot the estimated transferFunction with different methods

figure ; 
subplot(2,2,1) ; 
stem(h) ; 
title('Original Transfer function') ; 
subplot(2,2,2) ; 
stem(w_NLMS) ; 
title('Transfer function estimated with NLMS') ; 

subplot(2,2,3); 
stem(w_BLMS) ; 
title('Transfer function estimated with BLMS') ; 
subplot(2,2,4) ; 
stem(w_FDAF) ; 
title('Transfer function estimated with FDAF') ; 




% %% LMS Algo
% filterSize = 100; 
% x_padded = [zeros(filterSize ,1);x] ; 
% w = zeros(filterSize,1) ; 
% %i=filterSize ; 
% maxIter= signalLength -1 ; 
% %mu = 1 ; 
% error =0  ; 
% % i= filterSize ; 
% epsilon = 10^-10 ;  
% alpha = 0.99 ; 
% lambda = 0.1 ; 
% 
% i = 2*filterSize ; 
% N = filterSize; 
% while(i< signalLength - 2*filterSize -1 )
%              
%              %%non block LMS Algo
% %               x_n = x_padded(j+filterSize  :-1:j+1) ; 
% %               e = D(j) - x_n'*w  
% %               mu = 1/(x_n'*x_n)   ; 
% %               w = w + mu.* x_n*e ; 
% %               error = [error; e] ;
% %              
%              
%              
%              [x_matrix, autocorrX, mu ]=formX(filterSize, i, x_padded, N) ;
%             % block LMS Algo 
%             
%              E=eig(autocorrX)  
%              mu =alpha* 1/max(E) ;
%               modes = 1- 2*mu*E  
% %              mu = 1/(x_matrix*x_matrix' )  ; 
%              e = D(i-filterSize:i-filterSize + N -1)-x_matrix*w   ;   
%              error = [error ; e] ; 
%              w = w + mu*(1/N)*x_matrix'*e  ;  
%              
% 
%              
%              
% 
%              
%              if( (e'*e)/N< epsilon)
%                 
%               break ; 
%              end
%              
%      
%              
% i  = i + N ;   
% 
% end
% figure ; 
% plot(error) ; 
% title('error')  ; 
% 
% figure ; 
% stem(h, 'b')
% hold on ; 
% stem(w, 'r') ;
% title('System Identification by Adaptive LMS Algorithm')
% legend('Actual Filter Weights','Estimated Filter Weights',...
%        'Location','NorthEast')
% 
% figure ; 
% plot(D,'r') ; 
% hold on  ; 
% plot(filter(w,1,x),'b') ; 
% legend('Observed Signal d[n]','Filtred reference Signal y[n]',...
%        'Location','NorthEast')
% 
%    
% %% Filtring in the fourrier domain 
% 
% % dFFT = fft(D,2000) ; 
% % xFFT = fft(x,2000) ; 
% % figure ; 
% % plot(abs(fft(b,2000))) ; 
% % title('Magnitude of DFT with 2000 taps') ; 
% % 
% % w_fft = (dFFT./(xFFT)) ; 
% % figure ; 
% % plot(abs(w_fft)) ; 
% % title('Magnitude of the fourrier transform of the filter') ; 
% % 
% % figure  ; 
% % w_estimated = ifft(w_fft) ; 
% % stem(w_estimated) ; 
% % title ('estimated Filter')
% % 
% % figure  ; 
% % plot(D-filter(w,1,x)) ; 
% % title('error computed with estimated impulse response with fourier')
% % 
