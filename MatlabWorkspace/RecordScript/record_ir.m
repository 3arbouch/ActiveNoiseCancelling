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
   fs  =   96000;

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
   endTime = 5;
   sine= 0.5*sin(2*pi*f0/fs.*(1:endTime*fs)');
   sine = sine' ;

   %%% THIS ONE IS IMPORTANT! IT SET BASICALLY THE VOLUME!!!
   % Do NOT GO TOO HIGH!!! WITH THE MOTU, YOU CANNOT CHANGE THE VOLUME FOR
   % THE OUTPUTS VIA SYSTEM-VOLUME OR SOMETHING LIKE THAT!!!!
   output_sweep        =   0.50*sweep;
   %output_sweep        =   0.0.*sweep;

signalToPlay = sine ;
   
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
   pahandle    =   PsychPortAudio('Open', dev, 3,3,fs,[1 1],[],[],[],1);
   PsychPortAudio('FillBuffer', pahandle, signalToPlay);
   PsychPortAudio('GetAudioData', pahandle ,ceil(length(signalToPlay)/fs));
   PsychPortAudio('Start', pahandle);
   WaitSecs(size(signalToPlay,2)/fs);
   PsychPortAudio('Stop', pahandle);
    [signal]    =   PsychPortAudio('GetAudioData', pahandle ,[],[],[]);
    PsychPortAudio('Close');
    signal      =   signal';

% 
% IR = [];
%    for ind = 1:size(signalToPlay,2)
% % In the following, the deconvolution procedure
%     ind =  1;
%    len =   length(signalToPlay);
%    %IR  =   zeros(len,1);
% 
%    Y   = fft(signal(:,ind),(length(signal)+size(signalToPlay,2)-1));
%    H   = fft(signalToPlay',(length(signal)+size(signalToPlay,2)-1));
%    G   =   Y./(H);
% 
% % => impulse response
%    IR_tmp  =   ifft(G);
%    IR(:,ind) = IR_tmp(1:floor(sig_len.*fs));
%    
%    % plot Room implulse reponse 
%    
%    end
%    
%    %% PLots 
%    
%    IR = IR(:,ind);
%    IR = IR (1:fs);
%    figure ;
%    stem(IR);
   
   
   figure; 
   plot(signalToPlay)
   hold on 
   plot(signal(:,1))