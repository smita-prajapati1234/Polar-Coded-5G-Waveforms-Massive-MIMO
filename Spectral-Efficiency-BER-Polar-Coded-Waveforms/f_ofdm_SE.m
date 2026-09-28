clc;
clear all;
% close all;
s = rng(211);       % Set RNG state for repeatability

%Set parameters
numFFT = 1024;           % Number of FFT points
numRBs = 10;             % Number of resource blocks
rbSize = 20;             % Number of subcarriers per resource block
cpLen = 24;              % Cyclic prefix length in samples

bitsPerSubCarrier = 6;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM 10: 1024QAM
snrdB = 1:2:30;              % SNR in dB

toneOffset = 2.5;        % Tone offset or excess bandwidth (in subcarriers)
L = 513;                 % Filter length (=filterOrder+1), odd

%Filter Design
numDataCarriers = numRBs*rbSize;    % number of data subcarriers in sub-band
halfFilt = floor(L/2);
n = -halfFilt:halfFilt;

% Sinc function prototype filter
pb = sinc((numDataCarriers+2*toneOffset).*n./numFFT);

% Sinc truncation window
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;

% Normalized lowpass filter coefficients
fnum = (pb.*w)/sum(pb.*w);

% Filter impulse response
%  h = fvtool(fnum, 'Analysis', 'impulse', 'Fs', 15.36e6);

% Use dsp filter objects for filtering
filtTx = dsp.FIRFilter('Structure', 'Direct form symmetric', ...
    'Numerator', fnum);
filtRx = clone(filtTx); % Matched filter for the Rx
% Set up a figure for spectrum plot
hFig = figure(1);
% hFig = figure('Position', figposition([46 50 30 30]), 'MenuBar', 'none');
axis([-0.5 0.5 -200 50]);
hold on;
grid on
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['F-OFDM, ' num2str(numRBs) ' Resource blocks, '  ...
    num2str(rbSize) ' Subcarriers each'])

% Generate data symbols
bitsIn = randi([0 1], bitsPerSubCarrier*numDataCarriers, 1);

% QAM Symbol mapper
symbolsIn = qammod(bitsIn, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Pack data into an OFDM symbol
offset = (numFFT-numDataCarriers)/2; % for band center
symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
    zeros(numFFT-offset-numDataCarriers,1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));

% Prepend cyclic prefix
txSigOFDM = [ifftOut(end-cpLen+1:end); ifftOut];

% Filter, with zero-padding to flush tail. Get the transmit signal
txSigFOFDM = filtTx([txSigOFDM; zeros(L-1,1)]);

% Plot power spectral density (PSD)
[psd,f] = periodogram(txSigFOFDM, rectwin(length(txSigFOFDM)), ...
    numFFT*2, 1, 'centered');
plot(f,10*log10(psd),'b');

% Compute peak-to-average-power ratio (PAPR)
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
paprFOFDM = pm(txSigFOFDM);
disp(['Peak-to-Average-Power-Ratio for F-OFDM = ' num2str(paprFOFDM) ' dB']);


% Plot power spectral density (PSD) for OFDM signal
[psd,f] = periodogram(txSigOFDM, rectwin(length(txSigOFDM)), numFFT*2, ...
    1, 'centered');
hFig1 = figure('Position', figposition([46 15 30 30]));
plot(f,10*log10(psd));
grid on
axis([-0.5 0.5 -100 -20]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['OFDM, ' num2str(numRBs*rbSize) ' Subcarriers']);

% Compute peak-to-average-power ratio (PAPR)
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
paprOFDM = pm(txSigOFDM);
disp(['Peak-to-Average-Power-Ratio for OFDM = ' num2str(paprOFDM) ' dB']);

for loopx=1:length(snrdB)

% Add WGN
rxSig = awgn(txSigFOFDM, snrdB(loopx), 'measured');

% Receive matched filter
rxSigFilt = filtRx(rxSig);

% Account for filter delay
rxSigFiltSync = rxSigFilt(L:end);

% Remove cyclic prefix
rxSymbol = rxSigFiltSync(cpLen+1:end);

% Perform FFT
RxSymbols = fftshift(fft(rxSymbol));

% Select data subcarriers
dataRxSymbols = RxSymbols(offset+(1:numDataCarriers));

% Plot received symbols constellation
switch bitsPerSubCarrier
    case 2  % QPSK
        refConst = qammod((0:3).', 4, 'UnitAveragePower', true);
    case 4  % 16QAM
        refConst = qammod((0:15).', 16,'UnitAveragePower', true);
    case 6  % 64QAM
        refConst = qammod((0:63).', 64,'UnitAveragePower', true);
    case 8  % 256QAM
        refConst = qammod((0:255).', 256,'UnitAveragePower', true);
end
% constDiagRx = comm.ConstellationDiagram( ...
%     'ShowReferenceConstellation', true, ...
%     'ReferenceConstellation', refConst, ...
%     'Position', figposition([20 15 30 40]), ...
%     'EnableMeasurements', true, ...
%     'MeasurementInterval', length(dataRxSymbols), ...
%     'Title', 'F-OFDM Demodulated Symbols', ...
%     'Name', 'F-OFDM Reception', ...
%     'XLimits', [-1.5 1.5], 'YLimits', [-1.5 1.5]);
% constDiagRx(dataRxSymbols);

% Channel equalization is not necessary here as no channel is modeled

% BER computation
BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxBits = qamdemod(dataRxSymbols, 2^bitsPerSubCarrier, 'OutputType', 'bit', ...
    'UnitAveragePower', true);
ber = BER(bitsIn, rxBits);
ber1(loopx)=ber(1);

disp(['F-OFDM Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB) ' dB']);
end
T=[0:10:1200];
 for d=1:length(T)
     q=bitsPerSubCarrier*numFFT*T(d);
 r=(numFFT+cpLen+L)*T(d);
 SE_FOFDM=q./r;
     SE_FOFDM(d)=SE_FOFDM(1);
 end


 
 
% Restore RNG state
rng(s);

figure(100)
hold on ;
semilogy(snrdB,smooth(ber1),'g')
set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('BER')
title('BER vs SNR for F-OFDM Modulation')

figure(101)
hold on 
semilogy(T,SE_FOFDM,'--*g')
