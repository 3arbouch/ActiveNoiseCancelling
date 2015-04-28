% This script is used to run several experiments to estimate the headphone
% filter 
clear all
close all 


endTime=30;
f0=800;
frameSize=1024;

%% Generate signals
song=dsp.AudioFileReader('song.mp3');
noise=dsp.SignalSource(0.25*randn(endTime*song.SampleRate, 1), frameSize);
sine= dsp.SignalSource(0.5*sin(2*pi*f0/song.SampleRate.*(1:endTime*song.SampleRate)'), frameSize);



%% Setup M-AUDIO
H=dsp.AudioPlayer('DeviceName','ProFire 610','SampleRate', song.SampleRate, 'ChannelMappingSource','Property','ChannelMapping', 1);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'ProFire 610', ...
            'SampleRate', song.SampleRate, ...
            'OutputDataType','double',...
            'QueueDuration', 2,...
            'SamplesPerFrame', frameSize ,...
            'ChannelMappingSource','Property',...
            'ChannelMapping', [1 2]);


null_signal=zeros(frameSize,1);
% Contain the recording from the reference microhpne and the error
% microphone (inside the dummy head)
recordingsNoise=zeros(1,2);
noiseData= 0 ;
tic
while (toc<endTime)
    audio1 = step(song);
    audio2 = step(sine);
    audio3 = step(noise);
    
    % Produce white noise in the right speacker of the  headphone 
      step(H,  audio3(:,1));
    
    % Collect data from the microphones
     
     multichannelAudioFrame = step(AudioInput);
     recordingsNoise=[recordingsNoise;multichannelAudioFrame];
   

end


%% Release all the objects
release(song); % release the input file
release(sine);
release(noise);
release(AudioInput);
release(H);  % release the audio output device



%% Estimation of the filters S1(z) and s2(z)
noiseData = noise.Signal;
x = noiseData(1:length(recordingsNoise(:,1)));
d1 = recordingsNoise(:,1) ;
d2 = recordingsNoise(:,2) ;





% Compute the corrolation 
datacorr1 = xcorr(d1,x);
datacorr2 = xcorr(d2,x);

offset = 20 ;

[~, index1]= max(datacorr1) ;
[~, index2]= max(datacorr2) ;

% delay = song.SampleRate + 400 ;

startPoint1 = index1 -length(x) - offset;
startPoint2 = index2 -length(x) - offset;

d1 = d1(startPoint1:end) ; 
d2 = d2(startPoint2:end) ; 


% Compute the estimation of the filters
convergenceThreshold = 10^-6;
secondsToprocess = 20 ;
numberOfIterations = song.SampleRate*secondsToprocess ;
filterSize = 256 ;
blockSize= filterSize ;

numberOfIterationsF = numberOfIterations/filterSize -1; 
[errorFDAF1, MSerrorFDAF1, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, S1]=FDAFOSM(x,d1, filterSize,blockSize, numberOfIterationsF, convergenceThreshold)  ;

[errorFDAF2, MSerrorFDAF2, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, S2]=FDAFOSM(x,d2, filterSize,blockSize, numberOfIterationsF, convergenceThreshold)  ;

figure; 
semilogy(MSerrorFDAF1);
title('Learning curve of the first filter ');

figure; 
semilogy(MSerrorFDAF2);
title('Learning curve of the second filter ');



figure ;
stem(S1);

figure ;
stem(S2);


%% Cancelling sine wave 


H=dsp.AudioPlayer('DeviceName','ProFire 610','SampleRate', song.SampleRate, 'ChannelMappingSource','Property','ChannelMapping', 3);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'ProFire 610', ...
            'SampleRate', song.SampleRate, ...
            'OutputDataType','double',...
            'QueueDuration', 2,...
            'SamplesPerFrame', frameSize ,...
            'ChannelMappingSource','Property',...
            'ChannelMapping', [1 2]);


null_signal=zeros(frameSize,1);
% Contain the recording from the reference microhpne and the error
% microphone (inside the dummy head)
recordings=zeros(1,2);
sineData= 0 ;
musicData = 0;
tic
while (toc<endTime)
    audio1 = step(song);
    audio2 = step(sine);
    audio3 = step(noise);
    
    % Produce white noise in the right speacker of the  headphone 
      step(H,   audio2(:,1));
    
    % Collect data from the microphones
%      musicData = [musicData; audio1(:,1)] ; 
     multichannelAudioFrame = step(AudioInput);
     recordings=[recordings;multichannelAudioFrame];
   

end


%% Release all the objects
release(song); % release the input file
release(sine);
release(noise);
release(AudioInput);
release(H);  % release the audio output device


%% Processing: Removing the echo

% Remove the eho from the reference mic by supressing music from headphones
% refMicData = recordings(:,1) ;
% filtredMusicData= filter(S1, 1,musicData );
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
% figure;
% subplot(2, 2, 1);
% plot(filtredMusicData);
% subplot(2,2,2);
% plot(refMicData);
% subplot(2,2,3)
% plot(denoisedSignalRefMic);
% 
% figure;
% subplot(2, 2, 1);
% [STFT1,F1,T1] = spectrogram(filtredMusicData,1024,512,1024,song.SampleRate);
% imagesc(T1,F1,log10(abs(STFT1)))
% subplot(2,2,2);
% [STFT2,F2,T2] = spectrogram(refMicData,1024,512,1024,song.SampleRate);
% imagesc(T2,F2,log10(abs(STFT2)))
% subplot(2,2,3)
% [STFT3,F3,T3] = spectrogram(denoisedSignalRefMic,1024,512,1024,song.SampleRate);
% imagesc(T3,F3,log10(abs(STFT3)))
% 
% 
% 
% 
% %% Estimation of P(z)
xPrime = filter(S2,1,recordings(:,1) );
errorMicData = recordings(:,2);
 

 
 errorMicData = errorMicData(46500:end);
 xPrime = xPrime(46500:end);

 
 

% % %Align the filtred and recorded signal


% datacorr = xcorr(errorMicData,xPrime);
% [~, index]= max(datacorr) ;
% offset = 20 ;
% startPoint = index -length(errorMicData) - offset;
% errorMicData = errorMicData(startPoint:end);

 

%
 figure ; 
 plot(errorMicData);
 hold on 
 plot(xPrime);
 
 
%  
convergenceThreshold = 10^-5;
secondsToprocess = 25;
numberOfIterations = song.SampleRate*secondsToprocess ;
filterSize = 1024 ;
blockSize= filterSize ;

numberOfIterationsF = numberOfIterations/filterSize -1; 
[errorFDAF3, MSerrorFDAF3, timeOfConvergenceNLMSFDAF3,timeOfComputationFDAF3, W]=FDAFOSM(xPrime,errorMicData, filterSize,blockSize, numberOfIterationsF, convergenceThreshold)  ;


figure ;
semilogy(MSerrorFDAF3);
figure;
stem(W)

%% Cancelling process
denoisingSignal = filter(W,1,xPrime); 

denoisedSignal = errorMicData -  denoisingSignal ;

figure;
subplot(2, 2, 1);
[STFT1,F1,T1] = spectrogram(errorMicData,1024,512,1024,song.SampleRate);
imagesc(T1,F1,log10(abs(STFT1)))
subplot(2,2,2);
[STFT2,F2,T2] = spectrogram(denoisingSignal,1024,512,1024,song.SampleRate);
imagesc(T2,F2,log10(abs(STFT2)))
subplot(2,2,3);
[STFT2,F2,T2] = spectrogram(denoisedSignal,1024,512,1024,song.SampleRate);
imagesc(T2,F2,log10(abs(STFT2)))



