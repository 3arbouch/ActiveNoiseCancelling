close all ;
clear all ;
load('echoFilter');
load('headphoneFilter');
load('insideHeadhoneFilter');


endTime=30;
f0=800;
frameSize=1024;

%% Generate signals
song=dsp.AudioFileReader('song.mp3');
noise=dsp.SignalSource(0.25*randn(endTime*song.SampleRate, 1), frameSize);
sine= dsp.SignalSource(0.5*sin(2*pi*f0/song.SampleRate.*(1:endTime*song.SampleRate)'), frameSize);



%% Setup M-AUDIO
H=dsp.AudioPlayer('DeviceName','ProFire 610','SampleRate', song.SampleRate, 'ChannelMappingSource','Property','ChannelMapping', [1 3]);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'ProFire 610', ...
            'SampleRate', song.SampleRate, ...
            'OutputDataType','double',...
            'QueueDuration', 2,...
            'SamplesPerFrame', frameSize ,...
            'ChannelMappingSource','Property',...
            'ChannelMapping', [1 2]);


recordings=zeros(1,2);
sineData= 0 ;
musicData = 0;
tic
while (toc<endTime)
    audio1 = step(song);
    audio2 = step(sine);
    audio3 = step(noise);
    
    % Produce white noise in the right speacker of the  headphone 
      step(H,   [audio1(:,1) audio2(:,1) ]);
    
    % Collect data from the microphones
     musicData = [musicData; audio1(:,1)] ; 
     multichannelAudioFrame = step(AudioInput);
     recordings=[recordings;multichannelAudioFrame];
   

end



%% Release all the objects
release(song); % release the input file
release(sine);
release(noise);
release(AudioInput);
release(H);  % release the audio output device



%% Cancelling the music from the reference Mic: Echo
refMicData = recordings(:,1) ;
filtredMusicData= filter(S1, 1,musicData );

datacorr = xcorr(refMicData,filtredMusicData);
[~, index]= max(datacorr) ;

%Align the filtred and recorded signal
startPoint = index -length(filtredMusicData) +1 ;
refMicData = refMicData(startPoint:end);



filtredMusicData = filtredMusicData(1: length(refMicData));



denoisedSignalRefMic = refMicData-filtredMusicData;

figure;
subplot(2, 2, 1);
plot(filtredMusicData);
subplot(2,2,2);
plot(refMicData);
subplot(2,2,3)
plot(denoisedSignalRefMic);

figure;
subplot(3, 1, 1);
[STFT1,F1,T1] = spectrogram(refMicData,1024,512,1024,song.SampleRate);
imagesc(T1,F1,log10(abs(STFT1)))
title('Reference Signal')
subplot(3,1,2);
[STFT2,F2,T2] = spectrogram(filtredMusicData,1024,512,1024,song.SampleRate);
imagesc(T2,F2,log10(abs(STFT2)))
title('Echo canceller Signal')
subplot(3,1,3)
[STFT3,F3,T3] = spectrogram(denoisedSignalRefMic,1024,512,1024,song.SampleRate);
imagesc(T3,F3,log10(abs(STFT3)))
title('Echo cancelled Signal')


%% Cancelling 
denoisedSignalRefMic = filter(S2, 1, denoisedSignalRefMic);
denoisingSignal = filter (W, 1, denoisedSignalRefMic);
errorMicData = recordings(:,2);





datacorr = xcorr(errorMicData,denoisingSignal);
[~, index]= max(datacorr) ;

%Align the filtred and recorded signal
startPoint = index -length(errorMicData) +1 ;
errorMicData = errorMicData(startPoint:end);
errorMicData = errorMicData(1:length(denoisingSignal));



figure ; 
plot(errorMicData);
hold on 
plot(denoisingSignal);

cancelledSiganl = errorMicData - denoisingSignal ; 


%% Plots
figure;
subplot(3, 1, 1);
[STFT1,F1,T1] = spectrogram(errorMicData,1024,512,1024,song.SampleRate);
imagesc(T1,F1,log10(abs(STFT1)))
title('Listened Music with noise')
subplot(3,1,2);
[STFT2,F2,T2] = spectrogram(denoisingSignal,1024,512,1024,song.SampleRate);
imagesc(T2,F2,log10(abs(STFT2)))
title('cancelling Signal')
subplot(3,1,3)
[STFT3,F3,T3] = spectrogram(cancelledSiganl,1024,512,1024,song.SampleRate);
imagesc(T3,F3,log10(abs(STFT3)))
title('Cancelled Signal: Music witout noise')
