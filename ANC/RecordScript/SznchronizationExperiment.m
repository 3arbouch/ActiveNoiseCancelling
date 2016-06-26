close all 
% first you have to initialize the additional toolbox that allows for
% simultaneous realtime recording and playback


%% Initialization
   InitializePsychSound

%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   
   
   %%
%select the connected device
   dev =   3;

%set the sampling rate to be used
   fs  =   44100;

%set the length of the test-signal to be used
   sig_len = 20;


% The following code will produce a sine-sweep.

   fs1 =   1;
   fs2 =   fs./2;
   %fs2 = 5000;

   T   =   sig_len;
% Create the swept sine tone
   w1  =   2*pi*fs1;
   w2  =   2*pi*fs2;
   K   =   T*w1/log(w2/w1); %MATLAB notation: log(x)==ln(x)
   L   =   T/log(w2/w1);
   t   =   linspace(0,T-1/fs,fs*T);
   sweep   =   sin(K*(exp(t/L) - 1));
   
   f0 = 3000 ;
   f1 = 10000 ;
   endTime = 10;
   sine= 0.5*sin(2*pi*f0/fs.*(1:endTime*fs)');
   sine2= 0.5*sin(2*pi*f1/fs.*(1:endTime*fs)');

   sine = sine' ;
   sine2= sine2';
   %%% THIS ONE IS IMPORTANT! IT SET BASICALLY THE VOLUME!!!
   % Do NOT GO TOO HIGH!!! WITH THE MOTU, YOU CANNOT CHANGE THE VOLUME FOR
   % THE OUTPUTS VIA SYSTEM-VOLUME OR SOMETHING LIKE THAT!!!!
   output_sweep        =   0.50*sweep;
   %output_sweep        =   0.0.*sweep;
   
   noise=0.25*randn(endTime*fs,1) ;

    
[signal,fs]=audioread('song.mp3');
signal = signal(:,1);
signal = signal';
buffer(1) = PsychPortAudio('CreateBuffer', [], sine);
buffer(2) = PsychPortAudio('CreateBuffer', [], signal(1:fs/2));   



%%%% end of seep-generation %%%%

% the following code manages the playback and recording.
% For reference, please see the psychportaudio online-help

% from online help:
% pahandle = PsychPortAudio('Open' [, deviceid][, mode][,
% reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
%
% mode 3 == full duplex
% reqlatencyclass 3 == most aggressive for audio device
% AFTER FS: [1,5] defines: 1 output, 5 inputs
   
    silence = zeros(1,5*fs);
   
   i=1 ;
   bufferSize = 1 ;
   frameSizeSeconds = 0.5 ;
   frameSize = frameSizeSeconds*fs;
   pahandle    =   PsychPortAudio('Open', dev, 3,3,fs,[1 1],[],[],[],1);
   PsychPortAudio('FillBuffer', pahandle, buffer(1));
    PsychPortAudio('GetAudioData', pahandle ,bufferSize);
   PsychPortAudio('Start', pahandle, 0, [], 1);
%    status = PsychPortAudio('GetStatus', pahandle);
   
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
            WaitSecs(0.001);
        end
    end

    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(audiodata(1,:)) >= voicetrigger)); %#ok<MXFND>
        
    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
     offset = 10;
    recordedaudio = audiodata(:, idx-offset:end);
    
   
   
   
   
   completeSignal = recordedaudio ;
   indexRefillBuffer = 800;
   nextSampleStartIndex= 1;
   while(nextSampleStartIndex + fs/2 -1< length(sine))
       
       x = signal(nextSampleStartIndex:nextSampleStartIndex+fs/2 -1);
       PsychPortAudio('RefillBuffer', pahandle, [], x, indexRefillBuffer);
       d =   PsychPortAudio('GetAudioData', pahandle ,[],frameSizeSeconds,frameSizeSeconds,[]);
       completeSignal = [completeSignal,d];
       nextSampleStartIndex = nextSampleStartIndex +fs/2 ;
       indexRefillBuffer = indexRefillBuffer + fs/2 ;
   end

    PsychPortAudio('Close');




   %%
   figure; 
   plot(signal(1:length(sine)))
   hold on 
   plot(completeSignal)