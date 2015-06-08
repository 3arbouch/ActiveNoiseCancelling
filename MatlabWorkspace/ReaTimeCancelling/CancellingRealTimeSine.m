%% This script is used to cancel the noise coming from outside 
% There is three main steps:
%  1) Take the echo cancelled signal from the reference mic by filtering
%  the signal coming from the headphones by S2(z)
%  2) Filter the the echo canceled signal with P(z) and put play it in the
%  headphones

close all 
clear all 
%% Initialization
   InitializePsychSound
   
   data = load('experiment33Filters(1024)');
   h = firpm(1000, [0 0.001 0.005 1], [0 0 1 1]);

%% 
%search for audio playback and recording devices
   devices = PsychPortAudio('GetDevices' );
   
% time of processing
timeOfProcessing = 30 ;
 
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

   [music,~] = audioread('song.mp3');
   music = music(:,1);
   silence = zeros(length(music),1);
   music = [sine; silence];
   music = music' ;
 % Create the buffer with the white noise  
silence = zeros(1,2*fs);
NoiseTime = length(music)-length(silence) ;
f0 = 800;
sineNoise= 0.7*sin(2*pi*f0/fs.*(1:NoiseTime)');
sineNoise = sineNoise';
sineNoise = [silence, sineNoise];

 silence = zeros(1,length(music));
 buffer = PsychPortAudio('CreateBuffer', [], [silence;silence;sineNoise;music]);
  
 
 
 % set the buffer size for the recordings
   bufferSize = 30 ;
   
 % Set the ize of the frame of the recording in seconds
    %frameSizeSeconds = ((filterSize/16)/fs) ;
    

   frameSizeSeconds = 15;
 
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
    
   
   completeSignal = recordedaudio ;
   counter = length(recordedaudio);
   totalCancelledSignal  = 0 ;
   i =0 ;
  
     
       d =   PsychPortAudio('GetAudioData', pahandle ,[],frameSizeSeconds,frameSizeSeconds,[]);
       
      
          portion = music(1:length(recordedaudio)+length(d));
          cancellerSignal = filter(data.S1,1,portion);
          completeSignalRefMic = [recordedaudio(4,:),d(4,:)];
          completeSignalErrorMic= [recordedaudio(3,:),d(3,:)];
          cancelledSignal = completeSignalRefMic - cancellerSignal ;
          step = length(completeSignalRefMic);
     
       
        temporary = filter (data.S2,1,cancelledSignal);
        sineCancellerSignal2 = filter(data.P,1,cancelledSignal);
        sineCancellerSignal = filter(data.P,1,temporary);
%         completeSignalErrorMic = filter(h,1,completeSignalErrorMic);
%         completeSignalErrorMic = completeSignalErrorMic(floor(length(h)/2 + 1):end);
%         completeSignalErrorMic = [completeSignalErrorMic, zeros(1,length(sineCancellerSignal)- length(completeSignalErrorMic))];
        denoisedSignal = completeSignalErrorMic-sineCancellerSignal ;

       
     
       figure ;
       plot(completeSignalErrorMic);
       hold on 
       plot(sineCancellerSignal)

       
       figure;
    subplot(2, 1, 1);
    [STFT1,F1,T1] = spectrogram(completeSignalErrorMic,1024,512,1024,fs);
    imagesc(T1,F1,log10(abs(STFT1)))
    title('Listened Music with noise')
    subplot(2,1,2)
    [STFT3,F3,T3] = spectrogram(denoisedSignal,1024,512,1024,fs);
    imagesc(T3,F3,log10(abs(STFT3)))
    title('Cancelled Signal: Music witout noise')

       
       
       PsychPortAudio('Close');
       save('ExperimentWithMusicData30SineV2', 'portion','completeSignalRefMic','completeSignalErrorMic');

      
       
       
   

    PsychPortAudio('Close');
    save('canceller','sineCancellerSignal2') ;
    save ('recordedCanceller','sineCancellerSignal')
    save ('denoisedSignal','denoisedSignal')


    
    
    
  %% Processing offline  

% 
%  
%  
%  
%  
%  refMicData = completeSignal(1,:) ;
% filtredMusicData= filter(S1, 1,music );
% 
% datacorr = xcorr(refMicData,filtredMusicData);
% [~, index]= max(datacorr) ;
% 
% %Align the filtred and recorded signal
% startPoint = index -length(filtredMusicData) +1 ;
% refMicData = refMicData(startPoint:end);
% 
% 
% 
% filtredMusicData = filtredMusicData(1: length(refMicData));
% 
% 
% 
% denoisedSignalRefMic = refMicData-filtredMusicData;
% 
%  
% 
%  
%  %% plots
%  figure ;
%  plot(refMicData)
%  hold on 
%  plot(filtredMusicData);
%  
% 
%  
%     
%     
    