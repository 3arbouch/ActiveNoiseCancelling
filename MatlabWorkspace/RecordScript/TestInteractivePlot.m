figure
hAx(1) = subplot(211);
hLine(1) = line('XData',t, 'YData',nan(size(t)), 'Color','b', 'Parent',hAx(1));
xlabel('Time (s)'), ylabel('Amplitude')
hAx(2) = subplot(212);
hLine(2) = line('XData',f, 'YData',nan(size(f)), 'Color','b', 'Parent',hAx(2));
xlabel('Frequency (Hz)'), ylabel('Magnitude (dB)')
set(hAx, 'Box','on', 'XGrid','on', 'YGrid','on')