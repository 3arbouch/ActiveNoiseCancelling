x = noiseData ; 
d = recordings(:,1); 

x = x(44100*2:end) ; 
d = d(44100*2:end) ; 

signalLength = length(x) ; 
filterSize= 44100*2 ; 
blockSize = filterSize ; 
numberOfIterations = 44100*4 ;  

[errorNLMS, MSerrorNLMS, timeOfConvergenceNLMS, w_NLMS]=NLMS(x,d, filterSize, numberOfIterations)  ;


figure  ; 
semilogy(MSerrorNLMS) ; 

figure ; 
bar(timeOfConvergenceNLMS) ; 


figure  ; 
stem(w_NLMS) ; 