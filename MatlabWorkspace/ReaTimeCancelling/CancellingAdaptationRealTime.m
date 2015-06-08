% This script performs the adaptation of the filters and the cancelling
% process in real Time:

close all 
clear all 
%% Initialization
   InitializePsychSound
   
   % NumberOfIterations
   numberIter = 4000 ;
   
   data = load('experiment33Filters(1024)');
   
   % Truncate the filter and determine the frameSize 
  % Determine the time of filght in samples to fix the filter & frame size
   
   additionalDelay= 0; 
   P = data.P;
   [~, index]= max(abs(P));
   P = P(index-additionalDelay:end);
   
   sampleDelay = index ; 
   frameSizePerSample = floor((index -additionalDelay)/2)  
 % frameSizePerSample = 256 ;
   S2 = data.S2;
   if(length(S2) > frameSizePerSample)
   S2= S2(1:frameSizePerSample);
   end
      if(length(P) > frameSizePerSample)

   P= P(1:frameSizePerSample);
      end
 
 
 %Define the filter size as the frameSizeperSample
 filterSize = floor(frameSizePerSample) 

 
 % Define the parameters of the FDAF Algorithm 
  
  W0 = zeros(2*filterSize, 1) ;
  w = zeros(filterSize,1) ;
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
  
   past1 = zeros (1,filterSize) ;  
   past2 = zeros (1,filterSize) ;
   past3 = zeros (1,filterSize) ; 
   past4 = zeros (1,filterSize) ;
      
      
%% 
%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   
% time of processing
timeOfProcessing = 20 ;
 
%select the connected device (M-Audio)
   dev =   3;

%set the sampling rate to be used
   fs  =   16000;
   
   
   
%set the length of the played signal
    endTime = 10 ;
   
% Import the song 
sineTime = 0.5 ;
f0 = 800;
sine= 0.7*sin(2*pi*f0/fs.*(1:sineTime*fs)');

[music,~] = audioread('song.mp3');
music = music(:,1);
 
 silence = zeros(length(music),1);
 syncSignal = [sine ; silence] ; 
 
%    music = music' ;
 % Create the buffer with the white noise  

buffer = PsychPortAudio('CreateBuffer', [], syncSignal');
  

 
 % set the buffer size for the recordings
   bufferSize = 10 ;
   
 % Set the size of the frame of the recording in seconds: Corresponds to
 % half the distance of flight 
    frameSizeSeconds = ((frameSizePerSample)/fs) ;
    
 
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

    
    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(audiodata(1,:)) >= voicetrigger)); %#ok<MXFND>
        
    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
    offset = 0 ;
    if(idx > 10)
     offset = 10;
    end
   offset
%     recordedaudio = audiodata(:, idx-offset:end);
%      l2 = length(recordedaudio);
%     temp = ceil(index-frameSizePerSample);

%        completeSignal = recordedaudio(3,:) ;
%        i=1 ;
    
complete = zeros(1,numberIter*frameSizePerSample) ;
completeCanceller = zeros(1,numberIter*frameSizePerSample) ;

 i = 0 ;
 w1 = (1/(frameSizePerSample))*ones(1,frameSizePerSample) ;  
 % Get current position of the headphone buffer
    status = PsychPortAudio('GetStatus', pahandle);
    startIndexBuffer = status.ElapsedOutSamples   
   while(i< numberIter)
       
       
if (i == 0)
% First iteration: just save the two chunks of recorded data 

%% Get the first chunk of data 
[d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);

%construct the canceller Signal: Get the past and current chunk of recorded
%data from the reference Mic and filter them with the initial guess of the
%filter 
cancellerSignal =  filter(w1,1,[past4 ; d(2,:)]) ; 
% Extract the desired part of the filtred verion of the signal
cancellerSignal = cancellerSignal(filterSize+1:end); 

completeCanceller(i*filterSize+1:(i+1)*filterSize) = cancellerSignal ; 

% Put the canceller signal in the right position: startIndex + delay
% between mics
t = startIndexBuffer + sampleDelay ; 
PsychPortAudio('RefillBuffer', pahandle, [], cancellerSignal, t );

% Update past vectors 
past4 = d(2,:) ; 

%% Get the second chink of data
[d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);

% calculate the canceller signal 
cancellerSignal =  filter(w1,1,[past4 ; d(2,:)]) ; 
cancellerSignal = cancellerSignal(filterSize +1:end) ; 

completeCanceller((i+1)*filterSize+1:(i+2)*filterSize) = cancellerSignal ; 


% Put the canceller signal in the right position 
t = t +  frameSizePerSample ; 
PsychPortAudio('RefillBuffer', pahandle, [], cancellerSignal, t );

% update the past signals 

past3 = past4 ;
past4 = d(2,:) ; 

else 
%% First chunk of data ;      
[d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);

% Perform the adpatation 

temp = filter(S2,1,[past1;past2;past3]);
temp = temp(filterSize+1:end) ; 
% Pass the error signal and the two previous recorded blocks filtred by S2
[ W, w ] = FDAFOSMOneIteration( temp', d(1,:)', filterSize, P_k, g,k ,W0);

%Save the errorSignal 
complete(i*frameSizePerSample+1:(i+1)*frameSizePerSample) = d(1,:) ; 


% Construct the canceller Signal with the new estimated filter w
cancellerSignal =  filter(w1,1,[past4 ; d(2,:)]) ; 
cancellerSignal = cancellerSignal(filterSize +1:end) ; 

completeCanceller(i*filterSize+1:(i+1)*filterSize) = cancellerSignal ; 

% Put the canceller signal in the right position 
t = t + frameSizePerSample ; 
PsychPortAudio('RefillBuffer', pahandle, [], cancellerSignal, t );

% update the pasts 
past1 = past2 ;
past2 = past3 ; 
past3 = past4 ; 
past4 = d(2,:) ; 

% Update the filter 
W0 = W ; 

%% Second chiunk of data 

[d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);

% Perform the adpatation 
temp = filter(S2,1,[past1;past2;past3]);
temp = temp(filterSize+1:end) ; 
% Pass the error signal and the two previous recorded blocks filtred by S2
[ W, w ] = FDAFOSMOneIteration( temp', d(1,:)', filterSize, P_k, g,k ,W0);

complete((i+1)*frameSizePerSample+1:(i+2)*frameSizePerSample) = d(1,:) ;


% Construct the canceller Signal with the new estimated filter w
cancellerSignal =  filter(w1,1,[past4 ; d(2,:)]) ; 
cancellerSignal = cancellerSignal(filterSize +1:end) ; 

completeCanceller((i+1)*filterSize+1:(i+2)*filterSize) = cancellerSignal ; 

% Put the canceller signal in the right position 
t = t + frameSizePerSample ; 
PsychPortAudio('RefillBuffer', pahandle, [], cancellerSignal, t );

% update the pasts 
past1 = past2 ;
past2 = past3 ; 
past3 = past4 ; 
past4 = d(2,:) ; 

% Update the filter 
W0 = W ; 



    

    
end



        i = i+2 ;
           
   end
   
    PsychPortAudio('Close');

    
    
    
    
    
%% Plots

figure ;
plot(complete)
hold on 
plot(completeCanceller) ; 


figure ; 
stem(w);
    
    
