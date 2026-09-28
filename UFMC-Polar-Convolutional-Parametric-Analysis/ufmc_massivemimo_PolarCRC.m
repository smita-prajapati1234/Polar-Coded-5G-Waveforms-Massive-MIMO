 clc 
 clear all
%  close all
 s = rng(211);       % Set RNG state for repeatability
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
 prototypeFilter = chebwin(filterLen, slobeAtten);
% prototypeFilter = kaiser(filterLen, 2.5);
% prototypeFilter = hamming (filterLen);

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
frsz=m*n1
data11=randi([0 1],1160,1);
N3 = 116;
tempdata=[];
for clp=1:10

dataplr1=data11((clp-1)*N3+1:clp*N3);
dataplr = nrCRCEncode(dataplr1,'16');
% Find the power of 2 that is equal to or greater than N
n2 = ceil(log2(256));



% Polar encoding
encodedBits = nrPolarEncode(dataplr, 2^n2);
tempdata=[tempdata ;encodedBits(:)];
end


data1=[tempdata(:) ; zeros(frsz-rem(length(tempdata),frsz),1)];

for loopx=1:length(snrdBarr)


dataagg=[];
for dtl=1:ceil(length(data1)/frsz)
data=reshape(data1((dtl-1)*frsz+1:dtl*frsz),m,n1);
txSig=ufmc_mod(data);

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


% SETTING THE PARAMETERS FOR THE SIMULATION
fc=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/fc;        %Wavelength
d=l/2;        %Rx array spacing
N=64;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users
n=1:N;                          %Rx array number
n=transpose(n);                 %Row vector to column vector
sigma=0.1;

% RECEIVE SIGNAL MODEL (LINEAR)
H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);  %Channel matrix of size NxM
%  H=H/norm(H'*H);
H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);
w111=reshape([txSig;zeros(6,1)],M,35);
x=H*w111;                             %Receive vector of length N
% Add WGN
tempx=zeros(size(w111));


rxSig = awgn(x, snrdBarr(loopx),'Measured');

% LINEAR ARRAY PROCESSING - METHODS
% 1-MATCHED FILTER
 y=pinv(H'*H)*H'*rxSig;
% y=(H'*H+(2*sigma^2)*eye([M,M]))^(-1)*(H')*rxSig;
%  tempx=tempx+y;
% 
% y=tempx/1;


rxBits =ufmc_demod(y);

BER = comm.ErrorRate;

ber = BER(data(:), rxBits);

disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdBarr(loopx)) ' dB']);
ber_lp(dtl)=ber(1)

dataagg=[dataagg ;rxBits(:)];

end
tempdatad=[]
N1=2^n2;
for clp=1:10
% Polar decoding
decodedBits1 = nrPolarDecode((-100*(2*dataagg((clp-1)*N1+1:N1*clp)-1)), 132,256,8);
 decodedBits = nrCRCDecode(decodedBits1,'16');
tempdatad=[tempdatad; decodedBits(:)];
end

ber1(loopx)=mean(ber_lp)
[nerr berc(loopx)]=biterr(data11(:),tempdatad(:));
end
figure(100)

hold on
% semilogy(snrdBarr,smooth(ber1),'r')
semilogy(snrdBarr,smooth(berc),'Color','r')
set(gca,'Yscale','log')
xlabel('SNR')
ylabel('BER')
grid on 

