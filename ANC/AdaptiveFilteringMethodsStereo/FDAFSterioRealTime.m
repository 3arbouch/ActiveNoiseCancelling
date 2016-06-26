function [ distances, W1, W2, P_k1, P_k2] = FDAFSterioRealTime( signals, W1Previous, W2Previous, P_k1_previous, P_k2_previous,g,K)
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


 
x1 = signals(:,3) ;
x2 = signals(:,4);
d  = signals(1025:end,5) ; 

 
 % This is a fixed parameter, it's value was derived from experiments
 alpha= 0.3; 


  
 % This is a estimated parameter, it's value was derived from experiments:
 % Further research can be done on how to fix ths parameter 
 mu = 0.2 ; 
 






            X_k1= fft(x1) ;
            X_k2=fft(x2) ;

            
            Y_k =X_k1.*W1Previous + X_k2.*W2Previous ; 
            %Recover the last N samples correspending to linear convolution
            y_k = K*ifft(Y_k) ; 
            d_k = d; 
          
            % compute the error 
            e_k = d_k -y_k ;  
            % insert a zero block
            E_k = fft(K'*e_k) ; 
            
            P_k1 = (1-alpha)*P_k1_previous + alpha*abs((X_k1)).^2 ; 
            P_k_inv1 = 1./P_k1 ; 
            
            P_k2 = (1-alpha)*P_k2_previous + alpha*abs((X_k2)).^2 ; 
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
                
          
              
              W1 = W1Previous +  2*fft(gradient1) ;  
              W2 = W2Previous +  2*fft(gradient2) ; 
              
              
%               %Calculate the distances
%               w1 = ifft(W1) ; 
%               w1 = w1(1:filterSize) ;
%               w2 = ifft(W2) ; 
%               w2 = w2(1:filterSize) ;
%               
%                numberOfSamplesToIgnore = 50 ;
%                 [~, delay1] = max(abs(w1(numberOfSamplesToIgnore:filterSize-numberOfSamplesToIgnore)),[],2) ; 
%                 [~, delay2] = max(abs(w2(numberOfSamplesToIgnore:filterSize-numberOfSamplesToIgnore)),[],2) ; 
% 
%                 distance1 = delay1 +numberOfSamplesToIgnore ;
%                 distance2 = delay2 +numberOfSamplesToIgnore ;
%                 
                 distances =1 ; 


end










function [X_k] = formX( filterSize, index, data )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize-filterSize+1:1:index*filterSize+filterSize) ; 
X_k = fft(x_n) ;

 
end

