function [ error, MSerror, timeOfConvergence,timeOfComputation,  w ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,  numberOfIterations, convergenceThreshold )
% FDADOVM is the implementation of the FDAF algorithm based on the overlap
% save sectionning method.


% Returns:
% 1) The error computed at each sample
% 2) The mean sqaure error
% 3) Time of convergence: the time it takes to the MSE to be below the
%     convergence threshold
% 4) Time of computation: the time the algorithm takes for the whole
%   computaion
% 5) The estimated filter: the one that reaches the convergence threshold.
%    Otherwise, the filter estimated at the end of the computation


 
x = referenceSignal ; 
d = filtredSignal ; 

 
  W = zeros(2*filterSize, 1) ;

  P_k= ones(2*filterSize,1)  ; 
 g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
 
 % This is a fixed parameter, it's value was derived from experiments
 alpha= 0.3; 

 tic % Used to compute the computation time 
 
 error = 0 ; 
 MSerror = 0; 
 timeOfConvergence =0 ; 
  
 % This is a estimated parameter, it's value was derived from experiments:
 % Further research can be done on how to fix ths parameter 
 mu = 0.2 ; 
 
 i=1 ; 
 tstart = tic; % used to compute the convergence time
converged = 0 ;

 while (i<numberOfIterations)
            X_k=formX(filterSize, i, x) ;
            
            Y_k =X_k.*W ; 
            %Recover the last N samples correspending to linear convolution
            y_k = k*ifft(Y_k) ; 
            d_k = d(i*filterSize+1:(i+1)*filterSize) ; 
            
            % compute the error 
            e_k = d_k -y_k ;  
            % insert a zero block
            E_k = fft(k'*e_k) ; 
            
            P_k = (1-alpha)*P_k + alpha*abs((X_k)).^2 ; 
            P_k_inv = 1./P_k ; 
            %mu = alpha*(1/max(P_k)) ;
        
            % Set the step size: invertionally proportional to the power of
            % the signal
            mu_k = mu *P_k_inv ; 
            
            % compute the gradient 
            gradient = (conj(X_k).*E_k) ;
            % Append a zero block 
            GRADIENT = g*ifft(mu_k.*gradient);
            % Compute the update 
            W = W +  2*fft(GRADIENT) ; 
            
            error = [error ; e_k] ; 
            MSE = (error'*error)/length(error)  ;
             MSerror =  [MSerror ;repmat(MSE,filterSize,1) ] ; 
          %  MSerror =  [MSerror ;MSE ] ; 

               if(~converged && MSE < convergenceThreshold && i>2)
                 timeOfConvergence = toc(tstart) ;  
                  converged = 1;
                  W_convergence = W ;
               end
     
     
     i = i+ 1  ;
 end

timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ; 
        w = ifft(W) ; 
 else 
     w = ifft(W_convergence) ;
 end
           
w= w(1:filterSize) ; 



end










function [X_k] = formX( filterSize, index, data )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize-filterSize+1:1:index*filterSize+filterSize) ; 
X_k = fft((x_n)) ;

 
end

