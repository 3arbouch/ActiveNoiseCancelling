x= load('ListnedSoundWithNoise');
y = load('ListenedSoundCancelledNoise');

x= x.errorMicData;
y = y.cancelledSiganl ;
%sound(x(1:44100*10),44100);

sound(y(1:44100*10),44100);

