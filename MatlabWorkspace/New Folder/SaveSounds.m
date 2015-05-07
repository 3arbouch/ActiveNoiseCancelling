%% This script is used to run experiemnts and save them:
%  Signals needed to estimate S1(z) and S2(z): 
%    1) x[n]: reference signal (white noise) -> Generate white noise from 
%       speaker of the headphones
%    2) d1[n]:recorded signal at the reference mic (For S1(z)) -> Record
%       signal at reference mic
%    3) d2[n]: recorded signal at the error Mic  (For S2(z)) -> Record
%       signal at error Mic
% 
% Signals needed to estimate P(z):
%       Generate signal (white noise) at the 'Noise' loudspeacker    
%   1)  xPrime[n]: recorded signal at the reference mic : Record signal at the
%       reference mic
%   2)  dPrime[n]: recorded signal at the error mic: record signal at the error
%       microphone
% 

%% Initialization of the device and the signals: Generate ans Save signals for estimation of S1(z) and S2(z)


InitializePsychSound

%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   
   
 
%select the connected device (M-Audio)
   dev =   3;

%set the sampling rate to be used
   fs  =   44100;
   
%set the length of the played signal
    endTime = 10 ;
   
% Generate the white noise signal   
   noise=0.5*randn(endTime*fs,1) ;
   noise = noise';
   
 % Create the buffer with the white noise  
 
 buffer = PsychPortAudio('CreateBuffer', [], noise);
 
 % set the buffer size for the recordings
   bufferSize = 1 ;
   
 % Set the ize of the frame of the recording in seconds
   frameSizeSeconds = 0.5 ;

 
 % Open and start the device
 pahandle    =   PsychPortAudio('Open', dev, 3,3,fs,[1 2],[],[],[],1);
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
            level = max(abs(audiodata(2,:)));
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
    idx = min(find(abs(audiodata(2,:)) >= voicetrigger)); %#ok<MXFND>
        
    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
     offset = 5;
    recordedaudio = audiodata(:, idx-offset:end);
    
   completeSignal = recordedaudio ;
   counter = length(recordedaudio);
   tic
   while(counter< length(noise))
       
       d =   PsychPortAudio('GetAudioData', pahandle ,[],frameSizeSeconds,frameSizeSeconds,[]);
       completeSignal = [completeSignal,d];
       counter = counter + frameSizeSeconds*fs ; 
       
   end

    PsychPortAudio('Close');

    
    %% Plots
    figure  ;
    plot(noise);
    hold on ;
    plot(completeSignal(2,:));
    hold on 
    plot(completeSignal(1,:));
    
    
%% Save the recordings d1[n] and d2[n] and the refrecence signal x[n]
referenceMic = completeSignal(1,:);
errorMic = completeSignal(2,:);
save('referenceSignal','noise');
save('referenceMicSignal','referenceMic');
save('errorMicSignal','errorMic');


%% 
    

