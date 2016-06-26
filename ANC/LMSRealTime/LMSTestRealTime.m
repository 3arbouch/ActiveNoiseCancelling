close all
profile on 



convergenceThreshold = 10^-4;


%%  Generate signals of interest 
signalLength = 50000 ; 
originalFilterSize = 1024 ; 
[x,d,h] = generateSignals(signalLength, originalFilterSize,0 ) ; 






data = load('experiment50.mat');
filters = load ('experiment50Filters(1024)');



% filterSize =3* filterSize + 1 ; 
% d = d(delayPerSample+1:end);
filterSize = originalFilterSize ;
%% Parameters of the FDAFOSM algorithm
   W0 = zeros(2*filterSize, 1) ;
   w = zeros(filterSize,1) ;
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
  
  numberOfIterations = floor(length(x)/filterSize) - 1;
  completeMSerror = 0 ;
  for i =1:numberOfIterations -1
      % Take the past and the current block 
      x_n = x((i-1)*filterSize+1:1:(i+1)*filterSize) ;
      d_k = d(i*filterSize+1:(i+1)*filterSize) ;
    
      y_k = filter(w,1,x_n) ; 
     
      e_k = d_k - y_k(filterSize+1:end) ;
      
     
    
      [ W, w ] = FDAFOSMOneIteration( x_n, e_k, filterSize, P_k, g,k ,W0);
     
   
        completeMSerror = [completeMSerror,e_k'*e_k] ;
      W0 = W ;
  end


  % Parform the adapttaion with the hole signal
blockSize = filterSize ; 
convergenceThreshold = 10^-20 ;
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(x,d, filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;




  figure ;
  
  semilogy(completeMSerror);
  hold on ;
  semilogy(MSerrorFDAF)
  title ('Error with adaptation block by block')

 
figure ;
subplot(1,3,1);
stem(h) ; 
title ('Filter Impulse reponse') ;
subplot(1,3,2);
stem(w);
title ('Estimated Filter Impulse response with block adaptation') ;
subplot(1,3,3);
stem(w_FDAF);
title ('Estimated Filter Impulse response with hole signal adatation ') ;





% PLot the estimated transferFunction with different methods

