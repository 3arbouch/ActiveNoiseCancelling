function [ error , MSerror, timeOfConvergence, timeOfComputation , w ] = NLMS( referenceSignal, filtredSignal, filterSize, numberofIterations, convergenceThreshold )
% This function implement the NLMS
tic
x = referenceSignal ; 
d = filtredSignal ; 

i = filterSize ;
w = zeros(filterSize, 1) ; 
 NI =numberofIterations- i ;
% error = zeros(NI, 1) ;
 MSerror = zeros(NI- i, 1) ;
error = 0 ;
tstart = tic;
converged = 0 ;

while(i<numberofIterations)
               x_n = x(i:-1:i-filterSize + 1) ; 
               e = d(i) - x_n'*w  ;  
               mu = 1/(x_n'*x_n)   ; 
               w = w + mu.* x_n*e ; 
               j =i- filterSize +1;
              error = [error ; e] ; 
              MSerror(j) = (error'*error)/(j) ;
               if(~converged &&  MSerror(j) < convergenceThreshold && i > filterSize+1000)
                  timeOfConvergence = toc(tstart) ; 
                  converged = 1;
                  w_convergence = w ;
               end
              
              
               
  i = i+1  ;             
end
timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ;  
 else 
     w = w_convergence ;
 end

end

