function []=CancellationAdpatationOffline
% This script performs the adaptation of the filters and the cancelling
% process in real Time:

%%
% close all 
% clear all 
%    
% load('RealTimeCancellingExperiment')
%    data = load('experiment80Filters(4096, 16khz).mat');
% 
% 
% 
% referenceSignal = completeReferenceMic;
% filtredSignal = completeErrorMic;
% referenceSignal_prime= filter(data.S2,1,referenceSignal);
% filterSize = 4096 ; 
% blockSize = filterSize ; 
% 
% numberOfIterations = length(referenceSignal)/filterSize -5 ;
% convergenceThreshold = 10^-10 ; 
% 
% % numberOfIterations = 50 ;
% [ error, MSerror, timeOfConvergence,timeOfComputation,  P ] = FDAFOSM( referenceSignal_prime', filtredSignal', filterSize,blockSize,  numberOfIterations, convergenceThreshold );
% 
% figure ; 
% semilogy (MSerror)
% 
% figure ; 
% plot(completeErrorMic) ;
% 
% 
% figure ; 
% stem(P)
% 
% figure ; 
% cancellerSignal = conv(referenceSignal_prime, P) ; 
% cancellerSignal = cancellerSignal(1:length(filtredSignal)) ; 
% size(cancellerSignal(1:length(filtredSignal)));
% size(filtredSignal) 
%  plot(filtredSignal);  
%  hold on 
% plot( cancellerSignal);
%  
%  
%  cancelledSignal = completeErrorMic - cancellerSignal ; 
%%

   % Define the recording time in seconds 
   recordingTime = 20 ; 
   
   fs = 32000 ; 
   
 load('RealTimeCancellingExperiment2')

   
   filters = load('experiment90Filters(2048, 32khz)');
   
   % Truncate the filter and determine the frameSize 
  % Determine the time of filght in samples to fix the 
   
   P = filters.P;
   [~, index]= max(abs(P));
   

% Set the samples of flight as the maximum index of the offline estimated filter P(z)    
samplesTimeOfFlight = index ; 
sampleDelay = samplesTimeOfFlight/8  

% Take the frame Size as the nearest power of two to sample Delay in the
% lower band

power = (nextpow2(sampleDelay) -1) ; 
frameSize = 2^power    
   
% Set the filter Size
filterSize  = 2048 ;

% Calculate the number of frames that has to be recorded before the adaptation:
% We need to record filterSize samples to pass the recorded block to the
% FDAF Algorithm for adaptation
numberOfFrameBeforeAdaptation  = filterSize/frameSize  ; 



% Vector that will contain the output of the filtred blocks: the canceller
% signal 
lengthOutput = recordingTime*fs ;  

output = zeros(1, lengthOutput) ; 

% Truncate S2 if it's longer that the fixed filter Size
   S2 = filters.S2;
   size(S2)
%    S2 = S2(11:end) ; 
 
 
% Define vectors to store the past recorded signals in the reference mic
pastRefMic1 = zeros(1, filterSize) ; 
pastRefMic2 = zeros(1, filterSize) ; 
actualRefMic = zeros(1, filterSize) ;

% Define a vector to store the errorSignal recorded at the microphone 
actualErrorMic = zeros(1, filterSize) ;
 
% Define a vector that will be used for adaptation: Will contain the past vectors of the reference Signal


reference = zeros (1, 3*filterSize) ; 
 
 % Define the parameters of the FDAF Algorithm 
  
  W0 = zeros(2*filterSize, 1) ;
  %w = zeros(filterSize,1) ;
  w= P ; 
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
  
      
      

 
 % Set the size of the frame of the recording in seconds: PsychAudio getData method 
 % needs to set the frame in seconds
frameSizeSeconds = (frameSize/fs) ; 
    
 
 
 
% Vector to store the recorded error signal recorded at the error mic



 i = 0 ;
 j=0 ; 
  

   
    filters = zeros(filterSize, 1) ; 
    
 numberOfIterations = ((recordingTime-1)*fs)/frameSize -1 ; 
    
%  d = zeros(2,recordingTime*fs) ;
 MSError =  0 ; 
while(i<  numberOfIterations)
       
    % Acquire samplesToAddToReachFrameSize seconds obtain a frame Size 
d(2,:) = completeReferenceMic(i*frameSize+1:(i+1)*frameSize) ; 
d(1,:) = completeErrorMic(i*frameSize+1:(i+1)*frameSize) ; 
      

% save the recorded data from ref mic & error mic for Adaptation 
 actualRefMic(j*frameSize+1:(j+1)*frameSize) = d(2,:) ; 
 actualErrorMic(j*frameSize+1:(j+1)*frameSize) = d(1,:) ; 
 

 % Perform the convolution:  the filtering
 y_n = conv(w,d(2,:)) ; 
 
 %update the output vector 
[output,filtredBlock] =  updateOutput (output, y_n',i*frameSize+1 ) ; 

 %Take the corresconping block & put it in the right position of the
 %headphone buffer
 

j= j+1 ;  
i = i+1 ; 
t = t + frameSize ; 

    % This operation should tatke less that TimeOfFlight/2 !!!
    % Check if we collect filterSize samples to perform the adaptation
    if(j == numberOfFrameBeforeAdaptation)
        % construct the errorSignal 
      
%      tic
      % perform the filtering with S2 first ;      
      reference (1:filterSize) = pastRefMic2 ; 
      reference(filterSize+1:2*filterSize) = pastRefMic1 ; 
      reference(2*filterSize+1:end) = actualRefMic ; 
      temporary = filter(S2(1:filterSize/2), 1, reference) ; 
      temporary = temporary(filterSize+1:end) ;
      
       y_adpt = output((i-j)*frameSize+1:(i)*frameSize)  ;
       y_adpt = filter(S2(1:filterSize/2), 1,y_adpt);
       e_adpt = actualErrorMic - y_adpt ;
       MSError =[MSError,e_adpt] ; 
      
      [ W, w ] = FDAFOSMOneIteration( temporary', e_adpt', filterSize, P_k, g,k ,W0);
      filters = [filters, w] ; 
      pastRefMic2 = pastRefMic1 ; 
      pastRefMic1 = actualRefMic ; 
     
      W0= W ;  
      j = 0 ; 
%       toc
    end


    
end




    
    
    
    
    
%% Plots 



tfplot(w,fs,'w','EstimatedFilter real time');

tfplot(P,fs,'p','Filter estimated offline');

figure ; 
plot(completeErrorMic);
title('Signal recorded at error Mic')

figure ; 
plot(completeReferenceMic) ; 
title('Signal recorded at reference Mic')


figure ; 
plot(output)
title ('canceller Signal')

figure ; 
plot(MSError)
title ('errorSignal')



figure ; 
plot(completeErrorMic- output) ; 
hold on 
plot(output) ; 

dataCorr = xcorr(completeErrorMic,output); 
[value, index] = max(dataCorr) ; 

delay = index - length(output)
    
end

function [outputNew, filtredBlock] = updateOutput( output , y, index )
output(index: index+ length(y) -1 ) = output(index:index+ length(y) -1 ) + y ;  
 outputNew = output ; 
 filtredBlock = output(index: index+ length(y) -1 ); 
end

