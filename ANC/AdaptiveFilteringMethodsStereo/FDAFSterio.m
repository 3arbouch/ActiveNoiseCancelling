function [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = FDAFSterio( referenceSignal1, referenceSignal2, filtredSignal, filterSize,  convergenceThreshold )
% FDADOVM is the implementation of the sterio FDAF algorithm based on the overlap
% save sectionning method.
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
 
%Define the window
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

              
              W1 = W1 +  2*fft(window'.*gradient1) ;  
              W2 = W2 +  2*fft(window'.*gradient2) ; 
              
              
              
             error(i) = e_k'*e_k ; 
             MSE = (error*error')/i ;
             MSerror(i) = MSE ; 

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

