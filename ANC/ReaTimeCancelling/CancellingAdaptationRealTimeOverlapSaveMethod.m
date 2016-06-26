function []=CancellingAdaptationRealTimeOverlapSaveMethod
% This script performs the adaptation of the filters and the cancelling
% process in real Time:

close all 
clear all 
%% Initialization
   InitializePsychSound
   
   % Define the recording time in seconds 
   recordingTime = 60 ; 
   
   fs = 8000 ; 
   

   
   data = load('experiment100Filters(1024, 8khz)');
   
   % Truncate the filter and determine the frameSize 
  % Determine the time of filght in samples to fix the 
   
   P = data.P;
   [~, index]= max(abs(P));
   

% Set the samples of flight as the maximum index of the offline estimated filter P(z)    
samplesTimeOfFlight = index ; 
sampleDelay = samplesTimeOfFlight/8  

% Take the frame Size as the nearest power of two to sample Delay in the
% lower band

power = (nextpow2(sampleDelay) -1) ; 
frameSize = 2^power    
   
% Set the filter Size
filterSize  = 1024 ;

% Calculate the number of frames that has to be recorded before the adaptation:
% We need to record filterSize samples to pass the recorded block to the
% FDAF Algorithm for adaptation
numberOfFrameBeforeAdaptation  = filterSize/frameSize  ; 



% Vector that will contain the output of the filtred blocks: the canceller
% signal 
lengthOutput = recordingTime*fs ;  

output = zeros(1, lengthOutput) ; 

% Truncate S2 if it's longer that the fixed filter Size
   S2 = data.S2;
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
  
      
      
%% 
%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   

 
%select the connected device (M-Audio)
   dev =   3;

   
   
   
%set the length of the played signal
    endTime = 10 ;
   
% Import the song create the sync signal
sineTime = 0.5 ;
f0 = 800;
sine= 0.7*sin(2*pi*f0/fs.*(1:sineTime*fs)');

[music,~] = audioread('song.mp3');
music = music(:,1);
 
 silence = zeros(length(music),1);
 syncSignal = [sine ; silence] ; 
 

 % Create the buffer with the synchronization signal  

buffer = PsychPortAudio('CreateBuffer', [], syncSignal');
  

 
 % set the buffer size for the recordings
   bufferSize = 10 ;
   
 % Set the size of the frame of the recording in seconds: PsychAudio getData method 
 % needs to set the frame in seconds
frameSizeSeconds = (frameSize/fs) ; 
    
 
 % Open and start the device
 pahandle    =   PsychPortAudio('Open', dev, 3,3,fs,[1 2],[],[],[],1);
 PsychPortAudio('FillBuffer', pahandle, buffer);
 PsychPortAudio('GetAudioData', pahandle ,bufferSize);
 PsychPortAudio('Start', pahandle, 0, [], 1);
 
 % Synchronisation Step: Get the data once the sine is detected: Avoid
 % getting data of the M-Audio buffer
 
% Set the level and threshold  
   level = 0 ;
   voicetrigger = 10^-2; 
   
    while level < voicetrigger
        % Fetch current audiodata:
        [audiodata, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle);

        % Compute maximum signal amplitude in this chunk of data:
        if ~isempty(audiodata)
            level = max(abs(audiodata(1,:)));
        else
            level = 0;
        end
        
        % Below trigger-threshold?
        if level < voicetrigger
            % Wait for a millisecond before next scan:
            WaitSecs(0.00001);
        end
    end

     % Get current position of the headphone buffer
   
    
    
    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(audiodata(1,:)) >= voicetrigger)); %#ok<MXFND>
        
    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
    offset = 0 ;
%     if(idx > 10)
%      offset = 10;
%     end
   offset 
    recordedaudio = audiodata(:, idx-offset:end);
     l2 = length(recordedaudio)
     samplesToAddToReachFrameSize = frameSize - l2 ; 
     secondsToRecordToReachFramseSize = samplesToAddToReachFrameSize / fs ; 
     
     
%     temp = ceil(index-frameSizePerSample);

%        completeSignal = recordedaudio(3,:) ;
%        i=1 ;

% Vector to store the recorded error signal recorded at the error mic
completeErrorMic = zeros(1,recordingTime*fs) ;
completeReferenceMic = zeros(1,recordingTime*fs) ; 

% Vector to store the calculated canceller Signal
completeCanceller = zeros(1,recordingTime*fs) ;
Fblock = 0 ;

 i = 0 ;
 j=0 ; 
  
 % Get data as the filter has already output 
%  frameSizeBuffer = startIndexBuffer/fs
%   [d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeBuffer,frameSizeBuffer,[]);


%       t = startIndexBuffer + samplesTimeOfFlight  ;  
%      status = PsychPortAudio('GetStatus', pahandle);
%      startIndexBuffer = status.ElapsedOutSamples  

      t =   l2 ; 
    filters = zeros(filterSize, 1) ; 
    
 numberOfIterations = ((recordingTime-1)*fs)/frameSize -1 ; 
    
while(i<  numberOfIterations)
       
    % Acquire samplesToAddToReachFrameSize seconds obtain a frame Size 
%     if (i ==0 )
%         [r, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],secondsToRecordToReachFramseSize,secondsToRecordToReachFramseSize,[]);
%         d(2,:) = [recordedaudio(2,:),r(2,:)] ; 
%         d(1,:) = [recordedaudio(1,:),r(1,:)] ; 
%     else 
%         [d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);
%     end
%       
        [d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);

% save the recorded data from ref mic & error mic for Adaptation 
 actualRefMic(j*frameSize+1:(j+1)*frameSize) = d(2,:) ; 
 actualErrorMic(j*frameSize+1:(j+1)*frameSize) = d(1,:) ; 
 completeErrorMic(i*frameSize+1:(i+1)*frameSize) = d(1,:) ; 
 completeReferenceMic(i*frameSize+1:(i+1)*frameSize) = d(2,:) ; 

 % Perform the convolution:  the filtering
 y_n = conv(w,d(2,:)) ; 
 
 %update the output vector 
[output,filtredBlock] =  updateOutput (output,y_n',i*frameSize+1 ) ; 

 %Take the corresconping block & put it in the right position of the
 %headphone buffer
 
 PsychPortAudio('RefillBuffer', pahandle, [], filtredBlock, t );
 
%  [underflow, nextSampleStartIndex, nextSampleETASecs] = PsychPortAudio('FillBuffer', pahandle, filtredBlock , 1, t);


j= j+1 ;  
i = i+1 ; 
t = t + frameSize ; 

    % This operation should tatke less that TimeOfFlight/2 !!!
    % Check if we collect filterSize samples to perform the adaptation
    if(j == numberOfFrameBeforeAdaptation)
        % construct the errorSignal 
      
%       tic
      % perform the filtering with S2 first ;      
      reference (1:filterSize) = pastRefMic2 ; 
      reference(filterSize+1:2*filterSize) = pastRefMic1 ; 
      reference(2*filterSize+1:end) = actualRefMic ; 
      temporary = filter(S2, 1, reference) ; 
      temporary = temporary(filterSize+1:end) ;
      
%       temporary = [pastRefMic1, actualRefMic] ; 
      
      [ W, w ] = FDAFOSMOneIteration( temporary', actualErrorMic', filterSize, P_k, g,k ,W0);
%        filters = [filters, w] ; 
      pastRefMic2 = pastRefMic1 ; 
      pastRefMic1 = actualRefMic ; 
     
      W0= W ;  
      j = 0 ; 
%        toc
    end


    
end



    PsychPortAudio('Close');

    
    
    
    
    
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
plot(completeErrorMic) ;
hold on 
plot(output) ; 

dataCorr = xcorr(completeErrorMic,output); 
[value, index] = max(dataCorr) ; 

delay = index - length(output)
save('RealTimeCancellingExperimentReportV3') ; 
    
end

function [outputNew, filtredBlock] = updateOutput( output , y, index )
output(index: index+ length(y) -1 ) = output(index:index+ length(y) -1 ) - y ;  
 outputNew = output ; 
 filtredBlock = output(index: index+ length(y) -1 ); 
end

