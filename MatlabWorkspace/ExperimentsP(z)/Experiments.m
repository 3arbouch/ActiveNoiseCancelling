% This script is used to run several experiments to estimate the headphone
% filter 
clear all
close all 


endTime=12;
f0=7000;
frameSize=1024;

%% Generate signals
song=dsp.AudioFileReader('song.mp3');
noise= dsp.SignalSource(0.25*randn(endTime*song.SampleRate, 1), frameSize);
sine= dsp.SignalSource(0.5*sin(2*pi*f0/song.SampleRate.*(1:endTime*song.SampleRate)'), frameSize);



%% Setup M-AUDIO
H=dsp.AudioPlayer('DeviceName','ProFire 610','SampleRate', song.SampleRate, 'ChannelMappingSource','Property','ChannelMapping', 3);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'ProFire 610', ...
            'SampleRate', song.SampleRate, ...
            'OutputDataType','double',...
            'QueueDuration', 1,...
            'SamplesPerFrame', frameSize ,...
            'ChannelMappingSource','Property',...
            'ChannelMapping', [1 2]);


null_signal=zeros(frameSize,1);
recordings=zeros(1,2);
tic
while (toc<endTime)
    audio1 = step(song);
    audio2 = step(sine);
    audio3 = step(noise);
    
    % Produce white noise in the speackers 
    step(H,  audio3(:,1));
    
    
    multichannelAudioFrame = step(AudioInput);
    recordings=[recordings;multichannelAudioFrame];
   

end


%% Release all the objects
release(song); % release the input file
release(sine);
release(noise);
release(AudioInput);
release(H);  % release the audio output device



%% Play the recording on the macbook pro speakers
% recording=dsp.SignalSource(recordings(:,2),frameSize)
% H2=dsp.AudioPlayer('DeviceName','Built-in Output','SampleRate', song.SampleRate)
% 
% while ~isDone(recording)
%     audioFrame = step(recording);
%     step(H2,audioFrame);
%     
%     
% end
% 
% release(recording)
% release(H2)
%% Plot an STFT
% Reference mic
[STFT1,F1,T1] = spectrogram(recordings(:,1),1024,512,1024,song.SampleRate);
%Right error mic
% [STFT3,F3,T3] = spectrogram(recordings(:,3),1024,512,1024,song.SampleRate);
%Left error mic
[STFT2,F2,T2] = spectrogram(recordings(:,2),1024,512,1024,song.SampleRate);

figure ; 
imagesc(T1,F1,log10(abs(STFT1)))

figure ; 
imagesc(T2,F2,log10(abs(STFT2)))



