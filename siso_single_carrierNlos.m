clear all
close all
clc
% SISO single carrier single user downlink system

% System parameters
numTx = 1;
numRx = 1;
numUsers = 1;            
fc = 2e9;       
d0 = 1;
d = 10;
n = 2; % Path loss component
shadowing_genral = 8;

% Generate data stream
numBits = 1000;
bitStream = randi([0 1], 1, numBits);

% 16-QAM parameters
bitsPerSymbol = 4;
M = 16; % Modulation order

% Calculate number of symbols needed
numSymbols = ceil(numBits / bitsPerSymbol);

% Pad the bit stream with zeros if necessary
paddedBitStream = [bitStream, zeros(1, numSymbols * bitsPerSymbol - numBits)];

% Reshape bit stream into symbols (each row will be a symbol)
symbols = reshape(paddedBitStream, bitsPerSymbol, []).';

% Convert bits to decimal values
symbolIndices = bi2de(symbols, 'left-msb');

% QAM modulation (16-QAM)
modulatedSignal = qammod(symbolIndices, M, 'InputType', 'integer');

%transmitted signal power
P_tx=mean(abs(modulatedSignal).^2);
P_tx_dB=10*log10(P_tx);

% Path loss and shadowing
pathLoss_dB = 10 * n * log10(d / d0);
pathLoss_linear=10^(pathLoss_dB/10);
shadowing_dB = shadowing_genral * randn(1, numUsers);
shadowing_linear=10.^(shadowing_dB/10);
%total loss in db
totalLoss_dB=pathLoss_dB+shadowing_dB;
totalLoss_linear=10^(totalLoss_dB/10);

%apply path loss to the transmitted signal
s_pathAndShadowing=modulatedSignal/sqrt(totalLoss_linear);
%energy of the signal after path loss
pr_pathAndShadowing_db=10*log10(mean(abs(s_pathAndShadowing).^2));

% Rayleigh fading channel
h = (randn(size(s_pathAndShadowing)) + 1i * randn(size(s_pathAndShadowing))) / sqrt(2);

% Received signal including Rayleigh fading
receivedSignal = h .* s_pathAndShadowing;


% Add random AWGN without specifying SNR
noiseVariance = 0.5;
noise = sqrt(noiseVariance) * (randn(size(modulatedSignal)) + 1i * randn(size(modulatedSignal))) / sqrt(2);

%after signal has gone through path loss, shadowing and rayleigh plus noise
receivedSignal_afterNoise = receivedSignal + noise;

% Calculate the power of the received signal
receivedPower = 10*log10(mean(abs(receivedSignal).^2));

% Calculate the power of the noise
noisePower = 10*log10(mean(abs(noise).^2));

% Calculate the SNR in dB terms
snr_dB_calculated = receivedPower - noisePower;

% Display the results
disp(['Transmitted signal power(dB): ', num2str(P_tx_dB)]);
disp(['Path loss(dB): ', num2str(pathLoss_dB)]);
%disp(['Signal power(dB) after path loss: ', num2str(p_r_path_db)]);
disp(['Shadowing(dB): ', num2str(shadowing_dB)]);
disp(['Total loss(dB): ', num2str(totalLoss_dB)]);
disp(['Signal power after path loss and shadowing(dB): ', num2str(pr_pathAndShadowing_db)]);
disp(['Estimated noise power(dB): ', num2str(noisePower)]);
disp(['Received signal power(transmitted+shadowing+path loss+Rayleigh)(dB): ', num2str(receivedPower)]);
disp(['Calculated SNR (dB): ', num2str(snr_dB_calculated)]);
