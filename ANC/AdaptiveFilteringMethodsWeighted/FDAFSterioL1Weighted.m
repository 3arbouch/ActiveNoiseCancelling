function [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2,weights1, weights2, variances1,variances2] = FDAFSterioL1Weighted( referenceSignal1, referenceSignal2, filtredSignal, filterSize,  convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds )
% This function implements the L1 weighted Sterio version of FDAF algorithm.
% mu_k1, mu_k2 are the mean of the respective weighting function 
% interval1, interval2 represent the prior allwed error interval in
% samples: Usually set to filterSize
% Upper limit represents the upper limit of the weighting function 
% MemoryInSeconds represents how many previous weighting functions are
% taken into account fot the mean and variance calculation
% This function estimates two filters that represent the different paths
% between the referencesignal1, referencesignall2 and the reference signal

% Returns:
% 1) The error computed at each sample
% 2) The mean sqaure error
% 3) Time of convergence: the time it takes to the MSE to be below the
%     convergence threshold
% 4) Time of computation: the time the algorithm takes for the whole
%   computaion
% 5) The estimated filter: the one that reaches the convergence threshold.
%    Otherwise, the filter estimated at the end of the computation

fs =48000 ; 

 
x1 = referenceSignal1 ;
x2 = referenceSignal2 ;
d = filtredSignal ; 

 % Initialise parameters of the algotrithm 
 W1 = zeros(2*filterSize, 1) ;
 W2 = zeros(2*filterSize, 1) ;

 P_k1= ones(2*filterSize,1)  ; 
 P_k2= ones(2*filterSize,1)  ; 

 g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
 K = [zeros(filterSize, filterSize), eye(filterSize)] ; 
 
 % This is a fixed parameter, it's value was derived from experiments
 alpha= 0.3; 

 tic % Used to compute the computation time 
 

 timeOfConvergence =0 ; 
  
 % This is a estimated parameter, it's value was derived from experiments:
 % Further research can be done on how to fix ths parameter 
 mu = 0.2 ; 

 beta= 15 ; 


 
n= 1:2*filterSize ;
N=filterSize ; 
window= sin((n+0.5)*pi/N).^2 ;
 
 
 i=1 ; 
 tstart = tic; % used to compute the convergence time
 converged = 0 ;
 numberOfIterations = floor(length(referenceSignal1)/filterSize)  ; 
 filters1 = zeros(numberOfIterations,filterSize) ; 
 filters2 = zeros(numberOfIterations,filterSize) ; 
 error = zeros(1,numberOfIterations) ; 
 MSerror = zeros(1,numberOfIterations) ; 
 
 
 
weights1 = zeros(numberOfIterations,2*filterSize) ; 
weights2 = zeros(numberOfIterations,2*filterSize) ; 
 


variances1 = zeros(numberOfIterations,1) ; 
variances2 = zeros(numberOfIterations,1) ; 
 
 %% Form k as a gaussian distribution

sigma_k1 = interval1/2 ; 
sigma_k2 = interval2/2 ; 

x = [1:2*filterSize]' ;
norm1 = normpdf(x,mu_k1,sigma_k1);
norm2 = normpdf(x,mu_k2,sigma_k2);


memory = memoryInSeconds;
memorySamples = floor((memory*fs)/filterSize) ; 
mu_buffer1 = [((interval1/2).*randn(memorySamples,1) + mu_k1); zeros(numberOfIterations,1)];
mu_buffer2 = [((interval2/2).*randn(memorySamples,1) + mu_k2); zeros(numberOfIterations,1)];



k1 = (max(norm1) - norm1)*(upperLimit/max(norm1));
k2 = (max(norm2) - norm2)*(upperLimit/max(norm2));


samplesNoise = 100 ;
sigma = samplesNoise/2 ; 
noise = sigma*randn(memorySamples + 1,1);



 
 
 
% figure ; 
 while (i<numberOfIterations)
            X_k1=formX(filterSize, i, x1) ;
            
            X_k2=formX(filterSize, i, x2) ;

            
            Y_k =X_k1.*W1 + X_k2.*W2 ; 
            %Recover the last N samples correspending to linear convolution
            y_k = K*ifft(Y_k) ; 
            d_k = d(i*filterSize+1:(i+1)*filterSize) ; 
            
            % compute the error 
            e_k = d_k -y_k ;  
            % insert a zero block
            E_k = fft(K'*e_k) ; 
            
            P_k1 = (1-alpha)*P_k1 + alpha*abs((X_k1)).^2 ; 
            P_k_inv1 = 1./P_k1 ; 
            
            P_k2 = (1-alpha)*P_k2 + alpha*abs((X_k2)).^2 ; 
            P_k_inv2 = 1./P_k2 ; 
%             mu = alpha*(1/max(P_k)) ;
        
            % Set the step size: invertionally proportional to the power of
            % the signal
            mu_k1 = mu *P_k_inv1 ; 
            mu_k2 = mu *P_k_inv2 ; 

           
            
            % compute the gradient 
            GRADIENT1 = (conj(X_k1).*E_k) ;
            GRADIENT2 = (conj(X_k2).*E_k) ;


            % Append a zero block 
             gradient1 = g*ifft(mu_k1.*GRADIENT1);
             gradient2 = g*ifft(mu_k2.*GRADIENT2);
                
          
              % Save the computed filter
              w1 = ifft(W1) ; 
              filters1(i,:)= w1(1:filterSize) ;  
              w2 = ifft(W2) ; 
              filters2(i,:)= w2(1:filterSize) ; 

       
              
%              W1 = W1 +  2*fft(gradient1) + fft(- k1.*sign(w1)) ; 
%              W2 = W2 +  2*fft(gradient2) + fft(- k2.*sign(w2)) ; 
             
             W1 = W1 +  2*fft(window'.*(gradient1 - k1.*sign(w1))) ; 
             W2 = W2 +  2*fft(window'.*(gradient2 - k2.*sign(w2))) ; 


              
             error(i) = e_k'*e_k ; 
             MSE = (error*error')/i ;
             MSerror(i) = MSE ; 
             
             
             
             ignoreSamples = 70 ; 
               [~,mu_k1 ] = max(abs(w1(ignoreSamples:end-ignoreSamples)));
               [~,mu_k2 ] = max(abs(w2(ignoreSamples:end-ignoreSamples))) ;
               

               mu_k1 = mu_k1 + ignoreSamples ; 
               mu_k2 = mu_k2 +ignoreSamples ; 
                
                
               mu_buffer1(i+memorySamples) = mu_k1 ; 
               mu_buffer2(i+memorySamples) = mu_k2 ; 
       
               muVector1 = mu_buffer1(i:i+memorySamples) + noise  ; 
               muVector2 = mu_buffer2(i:i+memorySamples) + noise ; 
               
               mu_k1 = mean(muVector1)  ;
               sigma_k1 = std(muVector1) ; 
               
               
               mu_k2 = mean(muVector2)  ;
               sigma_k2 = std(muVector2) ;
        
               norm1 = normpdf(x,mu_k1,sigma_k1);
               norm2 = normpdf(x,mu_k2,sigma_k2);
               
               k1 = (max(norm1) - norm1)*(upperLimit/max(norm1));
               k2 = (max(norm2) - norm2)*(upperLimit/max(norm2));
               
               weights1(i,:) = k1 ; 
               weights2(i,:) = k2 ; 
             
               variances1(i) = sigma_k1 ; 
               variances2(i) = sigma_k2 ; 
             
        

               if(~converged && MSE < convergenceThreshold && i>2)
                 timeOfConvergence = toc(tstart) ;  
                  converged = 1;
                  W_convergence1 = W1 ;
                  W_convergence2 = W2 ;
               end
     
     
     i = i+ 1  ;
 end

timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ; 
        w1 = ifft(W1) ; 
        w2 = ifft(W2) ; 
 else 
     w1 = ifft(W_convergence1) ;
     w2 = ifft(W_convergence2) ;

 end
           
w1= w1(1:filterSize) ; 

w2 = w2(1:filterSize) ;


end










function [X_k] = formX( filterSize, index, data )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize-filterSize+1:1:index*filterSize+filterSize) ; 
X_k = fft(x_n) ;

 
end

