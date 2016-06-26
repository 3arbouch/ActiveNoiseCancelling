function [ error, MSerror, timeOfConvergence, w ] = FDAFOSMBLOCK( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations )
% FDADOVM is the implementation of the FDAF algorithm based on the overlap
% save sectionning

x = referenceSignal ; 
d = filtredSignal ; 

 
  W = zeros(filterSize + blockSize, 1) ;

  P_k= ones(filterSize + blockSize,1)  ; 
 g=[eye(filterSize),zeros(filterSize,blockSize); zeros(blockSize,filterSize), zeros(blockSize,blockSize)] ; 
 k = [zeros(blockSize, filterSize), eye(blockSize)] ; 
 alpha= 0.9; 
 mu = 1 ; 
 
 
 error = 0 ; 
 MSerror = 0; 
 epsilon = 10^-6 ; 
 timeOfConvergence =0 ; 

 i=filterSize/blockSize ; 
 tic
 while (i<numberOfIterations)
            X_k=formX(filterSize, i, x, blockSize) ;
            
            Y_k =X_k*W ; 
            %Recover the last N samples correspending to linear convolution
            y_k = k*ifft(Y_k) ; 
            d_k = d(i*blockSize+1:(i+1)*blockSize) ; 
            e_k = d_k -y_k ;  
            E_k = fft(k'*e_k) ; 
            
            P_k = (1-alpha)*P_k + alpha*abs(diag(X_k)).^2 ; 
            P_k_inv = 1./P_k ; 
            %mu = alpha*(1/max(P_k)) ;
            mu = 0.9 ; 
            mu_k = mu*diag(P_k_inv) ; 
            
            W = W +  2*fft(g*ifft(mu_k*ctranspose(X_k)*E_k)) ; 
            
            error = [error ; e_k] ; 
            MSE = (error'*error)/length(error)  
            MSerror =  [MSerror ;repmat(MSE,blockSize,1) ] ; 
            
               if(MSE < epsilon)
                  timeOfConvergence = toc ;  
                  break ; 
               end
     
     
     i = i+ 1  
 end

if(timeOfConvergence==0)
    timeOfConvergence = toc ; 
end
           

w = ifft(W) ; 
w= w(1:filterSize) ; 


end










function [X_k] = formX( filterSize, index, data, blockSize )
% This function construct two chunks of data (old and new) 

x_n = data(index*blockSize-filterSize+1:1:index*blockSize+blockSize) ; 
  X_k = diag(fft((x_n))) ; 
 
end

