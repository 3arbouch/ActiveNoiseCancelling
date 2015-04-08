function [ error, MSerror, timeOfConvergence, w ] = FDAFOSM( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations )
% FDADOVM is the implementation of the FDAF algorithm based on the overlap
% save sectionning

x = referenceSignal ; 
d = filtredSignal ; 

 
  W = zeros(2*filterSize, 1) ;

  P_k= ones(2*filterSize,1)  ; 
 g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
 alpha= 0.3; 

 
 
 error = 0 ; 
 MSerror = 0; 
 epsilon = 10^-10 ; 
 timeOfConvergence =0 ; 
 mu = 0.2 ; 
 i=1 ; 
 tic
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
            
%                if(MSE < epsilon)
%                   timeOfConvergence = toc ;  
%                   
%                end
     
     
     i = i+ 1  ;
 end

if(timeOfConvergence==0)
    timeOfConvergence = toc ; 
end
           

w = ifft(W) ; 
w= w(1:filterSize) ; 


end










function [X_k] = formX( filterSize, index, data, N )
% This function construct two chunks of data (old and new) 

x_n = data(index*filterSize-filterSize+1:1:index*filterSize+filterSize) ; 
X_k = fft((x_n)) ;

 
end

