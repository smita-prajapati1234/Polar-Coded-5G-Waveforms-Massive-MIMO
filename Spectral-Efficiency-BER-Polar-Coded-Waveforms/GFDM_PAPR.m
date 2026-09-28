clear all
close all
clc
%% Define simulation parameters
N = 1024; % Number of subcarriers
T = 1e-6; % Symbol time
Fs = N/T; % Sampling frequency
SNRdB = 10:2:30; % Range of SNR values to simulate
%% Generate random binary data
data = randi([0,1],1,N);
%% OFDM modulation
% Perform IFFT on data
x = ifft(data);
% Add cyclic prefix
cp_length = round(N/4);
x_cp = [x(end-cp_length+1:end), x];
% Serialize data
s = reshape(x_cp,1,[]);
%% GFDM modulation
% Generate random symbols
symbols = qammod(randi([0,3],1,N),4);
% Define prototype filter
h = rrc(T,Fs,'FilterSpanInSymbols',6,'RaisedCosineBeta',0.5);
% Apply pulse shaping
gfdm_signal = gfdm(symbols,h,T,Fs,'NumSubcarriers',N);




%% Simulate transmission over AWGN channel
for i = 1:length(SNRdB)
% Convert SNR from dB to linear scale
SNR = 10^(SNRdB(i)/10);
% Compute noise variance
noise_var = 1/SNR;
% Add noise to signals
rx_ofdm = awgn(s,SNRdB(i),'measured');
rx_gfdm = awgn(gfdm_signal,SNRdB(i),'measured');
% OFDM demodulation
rx_ofdm = reshape(rx_ofdm,length(s)/N,N+cp_length);
rx_ofdm(:,1:cp_length) = [];
rx_ofdm_data = fft(rx_ofdm,[],2);
rx_ofdm_data = reshape(rx_ofdm_data,1,[]);
% GFDM demodulation
rx_gfdm_data = gfdm_demod(rx_gfdm,h,T,Fs,'NumSubcarriers',N);
% Compute BER
[ber_ofdm(i),~] = biterr(data,rx_ofdm_data);
[ber_gfdm(i),~] = biterr(data,rx_gfdm_data);
% Compute PAPR
papr_ofdm(i) = 10*log10(max(abs(x_cp).^2)/mean(abs(x_cp).^2));
papr_gfdm(i) = 10*log10(max(abs(gfdm_signal).^2)/mean(abs(gfdm_signal).^2));
% Compute complexity
comp_ofdm(i) = N;
comp_gfdm(i) = N*6; % Assuming pulse shaping filter length of 6
end
%% Plot results
figure;
semilogy(SNRdB,ber_ofdm,'b--',SNRdB,ber_gfdm,'r-','LineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
legend('OFDM','GFDM');
figure;
plot(SNRdB,papr_ofdm,'b--',SNRdB,papr_gfdm,'r-','LineWidth',2);
xlabel('SNR (dB)');
ylabel('PAPR (dB)');
legend('OFDM','GFDM');
figure;
plot(SNRdB,comp_ofdm,'b--',SNRdB,comp_gfdm,'r-','LineWidth',2);
xlabel('SNR (dB)');
ylabel('Complexity');
legend('OFDM','GFDM');


function g = rrc(K, M, a)
% RRC - Return Root Raised Cosine filter (time domain)
t = linspace(-M/2, M/2, M*K+1);
t = t(1:end-1); t = t';
g = (sin(pi*t*(1-a))+4*a.*t.*cos(pi*t*(1+a)))./(pi.*t.*(1-(4*a*t).^2));
g(find(t==0)) = 1-a+4*a/pi;
g(find(abs(t) == 1/(4*a))) = a/sqrt(2)*((1+2/pi)*sin(pi/(4*a))+(1-2/pi)*cos(pi/(4*a)));

g = fftshift(g);
g = g / sqrt(sum(g.*g));
end