%% This script is used to cancel the noise coming from outside 
% There is three main steps:
%  1) Take the echo cancelled signal from the reference mic by filtering
%  the signal coming from the headphones by S2(z)
%  2) Filter the echo canceled signal with P(z) and put play it in the
%  headphones

close all 
clear all 
%% Initialization
   InitializePsychSound
   c = load('canceller.mat');
   dataS = load('denoisedSignal.mat');
   canceller = - c.sineCancellerSignal2 ;
   % create a 1024 taps moving average filter 
   
  h = firpm(1000, [0 0.001 0.005 1], [0 0 1 1]);
   
   % NumberOfIterations
   numberIter = 1000 ;
   
   data = load('experiment31Filters(1024)');
   
   % Truncate the filter and determine the frameSize 
   % Take an additional delay of 10 samples
   additionalDelay= 0; 
   P = data.P;
   [~, index]= max(abs(P));
   P = P(index-additionalDelay:end);
   
   frameSizePerSample = (index -additionalDelay)/2 ; 
   S2 = data.S2;
   if(length(S2) > frameSizePerSample)
   S2= S2(1:frameSizePerSample);
   end
      if(length(P) > frameSizePerSample)

   P= P(1:frameSizePerSample);
      end
   S2P = filter(S2,1,P);
%% 
%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   
% time of processing
timeOfProcessing = 20 ;
 
%select the connected device (M-Audio)
   dev =   3;

%set the sampling rate to be used
   fs  =   44100;
   
   
   
%set the length of the played signal
    endTime = 10 ;
   
% Import the song 
sineTime = 0.5 ;
f0 = 800;
sine= 0.5*sin(2*pi*f0/fs.*(1:sineTime*fs)');

sineTime2 = 2 ;
f1= 2000;
sine2 = 0.7*sin(2*pi*f1/fs.*(1:sineTime2*fs)');
   [music,~] = audioread('song.mp3');
   music = music(:,1);
   music2=music';
   silence = zeros(length(music),1);
   canceller(1:length(sine)) = sine ;
   music = canceller(1:fs*15);
%    music = music' ;
 % Create the buffer with the white noise  
silence = zeros(1,2*fs);
NoiseTime = length(music)-length(silence) ;
f0 = 800;
sineNoise= 0.7*sin(2*pi*f0/fs.*(1:NoiseTime)');
sineNoise = sineNoise';
sineNoise = [silence, sineNoise];
silence = zeros(1,length(sineNoise));
buffer = PsychPortAudio('CreateBuffer', [], [silence;silence;sineNoise;music]);
  

 
 % set the buffer size for the recordings
   bufferSize = 10 ;
   
 % Set the ize of the frame of the recording in seconds
    frameSizeSeconds = ((frameSizePerSample)/fs) ;
    
 
 % Open and start the device
 pahandle    =   PsychPortAudio('Open', dev, 3,3,fs,[4 4],[],[],[],1);
 PsychPortAudio('FillBuffer', pahandle, buffer);
 PsychPortAudio('GetAudioData', pahandle ,bufferSize);
 PsychPortAudio('Start', pahandle, 0, [], 1);
 
% Set the level and threshold  
   level = 0 ;
   voicetrigger = 10^-2; 
   
    while level < voicetrigger
        % Fetch current audiodata:
        [audiodata, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle);

        % Compute maximum signal amplitude in this chunk of data:
        if ~isempty(audiodata)
            level = max(abs(audiodata(3,:)));
        else
            level = 0;
        end
        
        % Below trigger-threshold?
        if level < voicetrigger
            % Wait for a millisecond before next scan:
            WaitSecs(0.00001);
        end
    end

%     status = PsychPortAudio('GetStatus', pahandle);
%     startIndexBuffer = status.ElapsedOutSamples
    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(audiodata(3,:)) >= voicetrigger)); %#ok<MXFND>
        
    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
    offset = 0 ;
    if(idx > 10)
     offset = 10;
    end
   offset
    recordedaudio = audiodata(:, idx-offset:end);
     l2 = length(recordedaudio);
    temp = ceil(index-frameSizePerSample);

       completeSignal = recordedaudio(3,:) ;
       i=1 ;
    
complete = zeros(1,numberIter*frameSizePerSample) ;
 
   while(i< numberIter)
       
       
 
[d, ~ ,~ ,~] = PsychPortAudio('GetAudioData', pahandle,[],frameSizeSeconds,frameSizeSeconds,[]);
completeSignal = [completeSignal , d(3,:)];
%     if(i==1)
%     tic
%      l1 = length(d) ;
%     silence = zeros(1,l1);
%     status = PsychPortAudio('GetStatus', pahandle);
%     startIndexBuffer = status.ElapsedOutSamples ;
%     % portion = music(l2+1:l2 + length(d) );
%     % cancellerSignal = filter(data.S1,1,portion);
%     % cancelledSignal = d(1,:) - cancellerSignal ;
%      portionToFilter = [zeros(1,l1), d(4,:)];
% %      temporary = filter (S2,1,portionToFilter);
%      cancellerSignal = filter(S2P,1,portionToFilter);
%      cancellerSignal = cancellerSignal(l1+1:end);
%      
% %        complete(l2+1:l2 + length(d) )= cancellerSignal ;
%        previousRecording = d(4,:);
% %     x =    [silence;silence;silence;cancellerSignal];
%      waitFor = frameSizeSeconds - toc 
% 
%      
%      
%      if(waitFor<0) 
%          waitFor = 0; 
%      end
%     t = startIndexBuffer  + temp - ceil(waitFor*fs) ;
% 
%     else 
%         t = t + l1;
%     % portion = music(l2+(i-1)*l1 + 1:l2 + i*l1);
%     % cancellerSignal = filter(data.S1,1,portion);
%     % cancelledSignal = d(1,:) - cancellerSignal ;
%       
%      portionToFilter = [previousRecording , d(4,:)];
% %      temporary = filter (S2,1,portionToFilter);
% 
%     cancellerSignal = filter(S2P,1,portionToFilter);
%     cancellerSignal = cancellerSignal(l1+1:end);
% %     complete(l2+(i-1)*l1 + 1:l2 + i*l1) = cancellerSignal ;
%      previousRecording = d(4,:);
% %     x =    [silence;silence;silence;cancellerSignal];
%     %  x= d(1,:)   ;

%     end


%    PsychPortAudio('RefillBuffer', pahandle, [], x, t );
        i = i+1 ;
           
   end
   
    PsychPortAudio('Close');

    
    
    %% Compare the signal recorded offline and the signal recorded online
    
    
    
    
    plot(completeSignal)
    hold on 
    plot(dataS.denoisedSignal);
     
    
    
