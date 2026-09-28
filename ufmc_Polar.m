 clc 
 clear all
%  close all
 s = rng(211);       % Set RNG state for repeatability
numFFT = 1024;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
  prototypeFilter = chebwin(filterLen, slobeAtten);
%  prototypeFilter = kaiser(filterLen, 2.5);
%  prototypeFilter = hamming (filterLen);

bitsPerSubCarrier = 8;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));
snrdBarr =0:2:30;              % SNR in dB


save initpara
m=bitsPerSubCarrier*subbandSize;
n1=numSubbands;
frsz=m*n1;
data11=randi([0 1],1320,1);
N3 = 132;
tempdata=[];
for clp=1:10

dataplr=data11((clp-1)*N3+1:clp*N3);
% Find the power of 2 that is equal to or greater than N
n2 = ceil(log2(N3));



% Polar encoding
encodedBits = nrPolarEncode(dataplr, 2^n2);
tempdata=[tempdata ;encodedBits(:)];
end


data1=[tempdata(:) ; zeros(frsz-rem(length(tempdata),frsz),1)];



dataagg=[];
for dtl=1:ceil(length(data1)/frsz)
data=reshape(data1((dtl-1)*frsz+1:dtl*frsz),m,n1);

for loopx=1:length(snrdBarr)
txSig=ufmc_mod(data);

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);

rxSig = awgn(txSig, snrdBarr(loopx),'Measured');

rxBits =ufmc_demod(rxSig);

BER = comm.ErrorRate;

ber = BER(data(:), rxBits)

disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdBarr(loopx)) ' dB']);
ber_lp(dtl)=ber(1)

dataagg=[dataagg ;rxBits(:)];

end
tempdatad=[]
N1=2^n2;
for clp=1:10
% Polar decoding
decodedBits = nrPolarDecode((-100*(2*dataagg((clp-1)*N1+1:N1*clp)-1)), 132,256,8);
tempdatad=[tempdatad; decodedBits(:)];
end

ber1(loopx)=mean(ber_lp)
[nerr berc(loopx)]=biterr(data11(:),tempdatad(:))
end

rng(s);
figure(100)

hold on
% semilogy(snrdBarr,smooth(ber1),'r')
semilogy(snrdBarr,smooth(berc),'r')
set(gca,'Yscale','log')
xlabel('SNR')
ylabel('BER')
grid on 

