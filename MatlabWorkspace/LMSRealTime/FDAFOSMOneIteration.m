function [W, w ] = FDAFOSMOneIteration( referenceSignal, errorSignal, filterSize, P_k, g,k ,W0)
% FDADOVM is the implementation of the FDAF algorithm based on the overlap
% save sectionning method 
% This script is used to performjust one iteration of the FDFOSM Algorithm
% implemented in the LMS folder: Will be used to update the filter in real
% time in each iteration.
% 

% Outputs:
% MSError : meanSquare error 
% w : the new calculated filter (time domain) 
% W : new calculated Filter in frequency domain 

% Inputs: 
% referenceSignal: signal captured at the reference microphone
% errorSignal = signal captured at the error Microphone 
% W0 = the old calculated filter
%  P,g,k : Utility matrices (check FDAFOSM function for more details)
% alpha : coeffecient that determine the pourcentage of energy to take between the differents slots of the signal 
%

x = referenceSignal ; 
e = errorSignal ; 

 
  alpha= 0.3; 
  mu = 0.05  ; 
     
            X_k = fft((x)) ;

            E_k = fft(k'*e) ; 
            
            P_k = (1-alpha)*P_k + alpha*abs((X_k)).^2 ; 
            P_k_inv = 1./P_k ; 
        
%             mu_k = mu*1./((0.001 +abs(X_k).^2));
%             mu_k = mu*sum(abs(X_k)/filterSize)^-1 *P_k_inv ; 
%             mu_k = mu*sum(abs(P_k).^2)^-1 *P_k_inv ; 
%                mu_k = ((sum(abs(P_k).^2)/(length(P_k)))^-1)*P_k_inv ; 

             mu_k = mu*P_k_inv ; 

            gradient = (conj(X_k).*E_k) ;
            GRADIENT = g*ifft(mu_k.*gradient);
            W = W0 +  2*fft(GRADIENT) ; 
            

w = ifft(W) ; 
w= w(1:filterSize) ; 



end











