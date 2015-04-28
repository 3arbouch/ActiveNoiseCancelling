function [ error, MSerror, timeOfConvergence,timeOfComputation,  w ] = FDAFCC(  referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold )
% This function implements the FDAF algorithm based on circular convlution 

tic
x = referenceSignal ; 
d = filtredSignal ; 

 
  W = zeros(filterSize, 1) ;

 P_k= ones(filterSize,1)  ; 
 alpha= 0.3; 

 
 
 error = 0 ; 
 MSerror = 0; 
 timeOfConvergence =0 ; 
 mu = 0.1 ; 
 i=1 ; 
 tstart = tic;
converged = 0 ;
 while (i<numberOfIterations)
            X_k=formX(filterSize, i, x, blockSize) ;
            D_k = fft(d(i*filterSize:(i+1)*filterSize-1)) ;
            Y_k =X_k.*W ; 
            E_k = D_k - Y_k ; 
            
            P_k = (1-alpha)*P_k + alpha*abs((X_k)).^2 ; 
            P_k_inv = 1./P_k ; 
            %mu = alpha*(1/max(P_k)) ;
        
            mu_k = mu *P_k_inv ; 
            GRADIENT = (conj(X_k).*E_k) ;
            W = W +  2*mu_k.*GRADIENT;
            
            error = [error ; ifft(E_k)] ; 
            MSE = (error'*error)/length(error)  ;
            MSerror =  [MSerror ;repmat(MSE,filterSize,1) ] ; 
            
               if(~converged && MSE < convergenceThreshold && i> 5)
                timeOfConvergence = toc(tstart) ;  
                  converged = 1;
                  W_convergence = W ;

               end
     
     
     i = i+ 1  ;
 end
 
timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ; 
      
 else 
      W = W_convergence ;
 end
           


w = ifft(W) ; 


end










function [X_k] = formX( filterSize, index, data, N )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize:1:(index+1)*filterSize - 1) ; 
X_k = fft((x_n)) ;



end

