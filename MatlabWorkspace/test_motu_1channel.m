clear all
close all

%% Test record from 1 channel with teh motu
Fs=44100;
nBits=24;
nChannels=4;
ID=2;

recorder = audiorecorder(Fs,nBits,nChannels,ID);


disp('Start speaking.')
recordblocking(recorder, 5);
disp('End of Recording.');


%%
play(recorder);


%% Get audio data
myRecording = getaudiodata(recorder);


%% Plot the waveform.
plot(myRecording);

%% Save Wav file

audiowrite('test_1mic.wav',myRecording,Fs,'BitsPerSample',24);

