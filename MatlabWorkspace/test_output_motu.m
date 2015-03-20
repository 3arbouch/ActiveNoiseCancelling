clear all
close all 
endTime=10;
f0=10000;
frameSize=1024;

%% Generate signals
song=dsp.AudioFileReader('song.mp3');
noise= dsp.SignalSource(0.25*randn(endTime*song.SampleRate, 1), frameSize);
sine= dsp.SignalSource(0.5*sin(2*pi*f0/song.SampleRate.*(1:endTime*song.SampleRate)'), frameSize);



%% Setup MOTU
H=dsp.AudioPlayer('DeviceName','MOTU Traveler','SampleRate', song.SampleRate, 'ChannelMapping', 1:3);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'MOTU Traveler', ...
            'SampleRate', song.SampleRate, ...
            'NumChannels', 3,...
            'OutputDataType','double',...
            'QueueDuration', 2,...
            'SamplesPerFrame', frameSize);


null_signal=zeros(frameSize,1);
recordings=zeros(1,3);
tic
while (toc<endTime)
    audio1 = step(song);
    audio2 = step(sine);
    audio3 = step(noise);
    step(H, [audio2(:,1), audio1(:,1), audio1(:,2)]);
    
    
    multichannelAudioFrame = step(AudioInput);
    recordings=[recordings;multichannelAudioFrame];
   

end
toc

%% Release all the objects
release(song); % release the input file
release(sine);
release(noise);
release(AudioInput);
release(H);  % release the audio output device



%% Play the recording on the macbook pro speakers
recording=dsp.SignalSource(recordings(:,3),frameSize)
H2=dsp.AudioPlayer('DeviceName','Built-in Output','SampleRate', song.SampleRate)

while ~isDone(recording)
    audioFrame = step(recording);
    step(H2,audioFrame);
    
    
end

release(recording)
release(H2)
%% Plot an STFT
% Reference mic
[STFT1,F1,T1] = spectrogram(recordings(:,1),1024,512,1024,song.SampleRate);
%Right error mic
[STFT3,F3,T3] = spectrogram(recordings(:,3),1024,512,1024,song.SampleRate);
%Left error mic
[STFT2,F2,T2] = spectrogram(recordings(:,2),1024,512,1024,song.SampleRate);

figure ; 
imagesc(T1,F1,log10(abs(STFT1)))

figure ; 
imagesc(T2,F2,log10(abs(STFT2)))

figure ; 
imagesc(T3,F3,log10(abs(STFT3)))

