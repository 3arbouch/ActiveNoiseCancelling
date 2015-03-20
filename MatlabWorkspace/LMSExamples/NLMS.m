function [ error , MSerror, timeOfConvergence , w ] = NLMS( referenceSignal, filtredSignal, filterSize, numberofIterations )
% This function implement the NLMS
x = referenceSignal ; 
d = filtredSignal ; 

i = filterSize ;
w = zeros(filterSize, 1) ; 
error = 0  ; 
MSerror = 0; 
epsilon = 10^-10 ; 
timeOfConvergence =0 ; 
tic
while(i<numberofIterations)
               x_n = x(i:-1:i-filterSize + 1) ; 
               e = d(i) - x_n'*w  ;  
               mu = 1/(x_n'*x_n)   ; 
               w = w + mu.* x_n*e ; 
              error = [error; e] ;
              MSE = (error'*error)/length(error)   ; 
              MSerror =  [MSerror ;MSE ] ; 
               if(MSE < epsilon)
                  timeOfConvergence = toc ;  
               end
              
              
               
  i = i+1               
end
 
if(timeOfConvergence==0)
    timeOfConvergence = toc ; 
end

end

