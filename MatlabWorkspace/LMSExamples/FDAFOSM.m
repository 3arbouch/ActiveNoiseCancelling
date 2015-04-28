function [ error, MSerror, timeOfConvergence,timeOfComputation,  w ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold )
% FDADOVM is the implementation of the FDAF algorithm based on the overlap
% save sectionning method 
 tic
x = referenceSignal ; 
d = filtredSignal ; 

 
  W = zeros(2*filterSize, 1) ;

  P_k= ones(2*filterSize,1)  ; 
 g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
 alpha= 0.3; 

 
 
 error = 0 ; 
 MSerror = 0; 
 timeOfConvergence =0 ; 
 mu = 0.2 ; 
 i=1 ; 
 tstart = tic;
converged = 0 ;
 while (i<numberOfIterations)
            X_k=formX(filterSize, i, x, blockSize) ;
            
            Y_k =X_k.*W ; 
            %Recover the last N samples correspending to linear convolution
            y_k = k*ifft(Y_k) ; 
            d_k = d(i*filterSize+1:(i+1)*filterSize) ; 
            e_k = d_k -y_k ;  
            E_k = fft(k'*e_k) ; 
            
            P_k = (1-alpha)*P_k + alpha*abs((X_k)).^2 ; 
            P_k_inv = 1./P_k ; 
            %mu = alpha*(1/max(P_k)) ;
        
            mu_k = mu *P_k_inv ; 
            gradient = (conj(X_k).*E_k) ;
            GRADIENT = g*ifft(mu_k.*gradient);
            W = W +  2*fft(GRADIENT) ; 
            
            error = [error ; e_k] ; 
            MSE = (error'*error)/length(error)  ;
             MSerror =  [MSerror ;repmat(MSE,filterSize,1) ] ; 
          %  MSerror =  [MSerror ;MSE ] ; 

               if(~converged && MSE < convergenceThreshold && i>5)
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










function [X_k] = formX( filterSize, index, data, N )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize-filterSize+1:1:index*filterSize+filterSize) ; 
X_k = fft((x_n)) ;

 
end

