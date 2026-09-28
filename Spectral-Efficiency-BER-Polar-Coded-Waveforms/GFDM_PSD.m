close all
clear all
clc
trans = transmitter();
trans.mapper('QAM');
% OFDM_signal = trans.modulate('OFDM',0,64);
GFDM_signal = trans.modulate('GFDM',0,16);
f = linspace(-trans.no_of_total_subcarriers/2, trans.no_of_total_subcarriers/2, 2*length(GFDM_signal)+1); f = f(1:end-1)';
%  plot(f, mag2db(flip(abs(fft(OFDM_signal, 2*length(OFDM_signal)))))/2, 'b');
%  plot(f, mag2db((abs(fft(OFDM_signal, 2*length(OFDM_signal)))))/2, 'r');
 hold on;
 plot(f, mag2db((abs(fft(GFDM_signal, 2*length(GFDM_signal)))))/2, 'r');
 hold on;
 
 GFDM_signal = trans.modulate('GFDM',0,16);
  plot(f, mag2db((abs(fft(GFDM_signal, 2*length(GFDM_signal)))))/2, 'g');
GFDM_signal = trans.modulate('GFDM',0,32);
grid on;
plot(f, mag2db((abs(fft(GFDM_signal, 2*length(GFDM_signal)))))/2, 'b');
GFDM_signal = trans.modulate('GFDM',0,64);
  plot(f, mag2db((abs(fft(GFDM_signal, 2*length(GFDM_signal)))))/2, 'c');
 GFDM_signal = trans.modulate('GFDM',0,128);
  plot(f, mag2db((abs(fft(GFDM_signal, 2*length(GFDM_signal)))))/2, 'y');
 ylim([-5, 35]);
xlim([-130, 130])
xlabel('Frequency'); ylabel('Power Spectral Density [dB]');
legend({'GFDM - 32 blocks', 'GFDM - 16 blocks', 'GFDM - 8 blocks','GFDM - 4 blocks','GFDM - 2 blocks'});

trans.no_of_blocks