clear all
close all
clc
%% Definition of system parameters
numTx = 2;        % Number of transmit antennas
numRx = 2;        % Number of receive antennas
numUsers = 500;     % Number of users
fc = 5e9;         % Carrier frequency
lambda = physconst('LightSpeed') / fc;
distances = randi([100, 300], numUsers, 1);  % Generate random distances between 100 and 300 meters
meanShadowing_dB = 0;      % Mean shadowing in dB
stdShadowing_dB = 4;       % Standard deviation of shadowing in dB
numRealizations = 1000;    % Number of Monte Carlo realizations
%% Desired Transmit Power
desiredPower_dBm = 33; % Desired transmit power in dBm
desiredPower_Watts = 10^((desiredPower_dBm - 30) / 10); % Convert dBm to Watts
%% Bit Generation
numBits=1000;
bit_stream_tx = randi([0 1], numBits,numTx);
M=2; 
modulated_tx =modulation(bit_stream_tx,'psk',M);
modulated_power=bandpower(modulated_tx);
for tx=1:numTx
    xt1(:,tx)=modulated_tx(:,tx)/sqrt(modulated_power(:,tx));
end
% signal
tx = phased.Transmitter('PeakPower',1,'Gain',3);
xt = tx(xt1);
numbSymbols=length(xt(:,1));
%% Show transmit power
for tx = 1:numTx
    disp(['Transmit power for Tx Antenna ', num2str(tx), ': ', num2str(pow2db(bandpower(xt(:,tx)))+30), ' dBm']);
end
%% Path Loss Initialization
pathLoss_dB = fspl(distances, lambda);
tx_pathloss=zeros(numbSymbols,numTx,numUsers);
for userIdx = 1:numUsers
    pathLossFactor = 10^(-pathLoss_dB(userIdx) / 10);
    for tx = 1:numTx
        tx_pathloss(:, tx, userIdx) = xt(:,tx) * sqrt(pathLossFactor);
    end
disp(['Power after path loss for user ',num2str(userIdx),': ',num2str(pow2db(bandpower(tx_pathloss(:,1,userIdx)))+30),' dBm']);
end
%% Shadowing
%Generate shadowing variations for each user
shadowing_dB = stdShadowing_dB * randn(numRealizations, numUsers, numTx) + meanShadowing_dB;
tx_shad_power = zeros(numUsers, numTx);
tx_shadowing = zeros(numbSymbols, numTx, numUsers);
for userIdx = 1:numUsers
    for tx = 1:numTx
        total_power = 0;
        for realization = 1:numRealizations
            shadowingFactor = 10^(-shadowing_dB(realization, userIdx, tx) / 20);
            shadowing_effect = tx_pathloss(:, tx, userIdx) * sqrt(shadowingFactor);
            % Accumulate the signals after shadowing
            tx_shadowing(:, tx, userIdx) = tx_shadowing(:, tx, userIdx) + shadowing_effect;
        end
        % Average signal after shadowing
        tx_shadowing(:, tx, userIdx) = tx_shadowing(:, tx, userIdx) / numRealizations;

        % Calculate the power of the signal after shadowing
        tx_shad_power(userIdx,tx) = bandpower(tx_shadowing(:,tx,userIdx));
        disp(['Power after path loss and shadowing for Tx Antenna ', num2str(tx), ' for User ', num2str(userIdx), ': ', num2str(pow2db(tx_shad_power(userIdx,tx))+30), ' dBm']); 
    end
end
%% Reyleigh affect
rayleigh_fading = (randn(numRx, numTx) + 1j * randn(numRx, numTx)) / sqrt(2);

%% Apply Rayleigh Fading
received_signal = zeros(numbSymbols, numTx, numUsers);
for userIdx = 1:numUsers
    % Apply Rayleigh fading
    faded_signal = zeros(numbSymbols, numTx);
    for rx = 1:numRx
        for tx = 1:numTx
            faded_signal(:, rx) = faded_signal(:, rx) + rayleigh_fading(rx, tx) * tx_shadowing(:, tx, userIdx);
        end
    end
    received_signal(:, :, userIdx) = faded_signal;

    % Initialize a variable to store the total power received for the user
    total_power_received_Watts = 0;
    % Calculate and display the power at the receiver in dBm
    for rx = 1:numRx
        power_received_Watts = mean(abs(received_signal(rx, :, userIdx)).^2);
        total_power_received_Watts = total_power_received_Watts + power_received_Watts;
        power_received_dBm = 10 * log10(power_received_Watts) + 30;
        disp(['Received Power at Rx Antenna ', num2str(rx), ' for User ', num2str(userIdx), ': ', num2str(power_received_dBm), ' dBm']);
    end

    % Calculate the average power received for the user
    average_power_received_Watts = total_power_received_Watts / numRx;
    average_power_received_dBm = 10 * log10(average_power_received_Watts) + 30;
    disp(['Average Received Power for User ', num2str(userIdx), ': ', num2str(average_power_received_dBm), ' dBm']);
end

%% Add AWGN
noise_variance = 0.5;
awgn_noise = sqrt(noise_variance / 2) * (randn(numbSymbols, numRx, numUsers) + 1j * randn(numbSymbols, numRx, numUsers));

%% Calculate Noise Power in dBm
noise_power_Watts = noise_variance;
noise_power_dBm = 10 * log10(noise_power_Watts) + 30;
%% Calculate SNR
snr_dB = zeros(numRx, numUsers);
received_with_noise=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for rx = 1:numRx
        % Calculate noise power
        noise_power_Watts = noise_variance;
        % Add AWGN to the received signal
        received_with_noise(:,rx,userIdx) = received_signal(:, rx, userIdx) + awgn_noise(:, rx, userIdx);      
        % Calculate received signal power
        signal_power_Watts = bandpower(received_signal(:,rx,userIdx));
        
        % Calculate SNR at each receiver
        snr_dB(rx, userIdx) = 10 * log10(signal_power_Watts / noise_power_Watts);
    end
end
%% Average SNR Calculation
average_snr_dB = mean(snr_dB, 1);
for userIdx = 1:numUsers
    disp(['Average SNR for User ', num2str(userIdx), ': ', num2str(average_snr_dB(userIdx)), ' dB']);
end
%% Calculate BER for BPSK
% Define the Q-function for BER calculation
Q = @(x) 0.5 * erfc(x / sqrt(2));

% Initialize BER and Throughput matrices
ber = zeros(numRx, numUsers);
throughput = zeros(numUsers, 1);

% Calculate BER for each RX
for userIdx = 1:numUsers
    for rx = 1:numRx
        % Assuming BPSK modulation
        ber(rx, userIdx) = Q(sqrt(10^(snr_dB(rx, userIdx) / 10)));
    end
end
% Average BER Calculation
average_ber = mean(ber, 1);
% data rate in bpsk=1
data_rate=1;
%% Calculate Throughput for each user
throughput = data_rate * (1 - average_ber);

% Display average BER and throughput for each user
for userIdx = 1:numUsers
    disp(['Average BER for User ', num2str(userIdx), ': ', num2str(average_ber(userIdx))]);
    disp(['Throughput for User ', num2str(userIdx), ': ', num2str(throughput(userIdx)), ' bits per second']);
end

%% Graphs
% SNR vs. Distance
figure;
scatter(distances, average_snr_dB, 'filled');
xlabel('Distance (meters)');
ylabel('Average SNR (dB)');
title('Average SNR vs. Distance');
grid on;

% Throughput vs. SNR
figure;
scatter(average_snr_dB, throughput, 'filled');
xlabel('Average SNR (dB)');
ylabel('Throughput (bits per second)');
title('Throughput vs. Average SNR');
grid on;
%% RIS setup
% Setup surface
Nr = 10;
Nc = 20;
dr = 0.5*lambda;
dc = 0.5*lambda;

% construct surface
ris = helperRISSurface('Size',[Nr Nc],'ElementSpacing',[dr dc],...
    'ReflectorElement',phased.IsotropicAntennaElement,'OperatingFrequency',fc);

dv = -1;
du = sqrt((distances.^2)-dv^2);
pos_ue = [du.'; dv*ones(1,numUsers); zeros(1,numUsers)];

% scene
dbr1 = 100; %for ris1
pos_ap = [0;0;0];
pos_ris1 = [dbr1;0;0]; 
v = zeros(3,1);

dbr2=dv;
pos_ris2=[0;dbr2;0];

% compute the range and angle of the RIS from the base station and the UE
[r_ap_ris,ang_ap_ris] = rangeangle(pos_ap,pos_ris1);
[r_ue_ris,ang_ue_ris] = rangeangle(pos_ue,pos_ris1);

fs = 10e6;
c = physconst('lightspeed');

% channel
chanAPToRIS = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);
chanRISToUE = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);

N0dB=-4;
rcoeff_ris = ones(Nr*Nc,1);
x_ris_in1=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_in1(:,i,userIdx) = chanAPToRIS(xt(:,i),pos_ap,pos_ris1,v,v);
    end
end
x_ris_in2=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_in2(:,i,userIdx) = chanAPToRIS(xt(:,i),pos_ap,pos_ris2,v,v);
    end
end

x_ris_out1=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_out1(:,i,userIdx) = ris(x_ris_in1(:,i,userIdx),ang_ap_ris,ang_ue_ris(:,userIdx),rcoeff_ris);
    end
end
x_ris_out2=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_out2(:,i,userIdx) = ris(x_ris_in2(:,i,userIdx),ang_ap_ris,ang_ue_ris(:,userIdx),rcoeff_ris);
    end
end

ynlos=received_signal;

ylos1ris=zeros(numbSymbols,numRx,numUsers);
SNRlos1ris=zeros(numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        ylos1ris(:,i,userIdx) = chanRISToUE(x_ris_out1(:,i,userIdx),pos_ris1,pos_ue(:,userIdx),v,v) + ynlos(:,i,userIdx);
        SNRlos1ris(i,userIdx) = pow2db(bandpower(ylos1ris(:,i,userIdx))) - N0dB;
    end
end

ylos2ris=zeros(numbSymbols,numRx,numUsers);
SNRlos2ris=zeros(numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        ylos2ris(:,i,userIdx) = chanRISToUE(x_ris_out1(:,i,userIdx),pos_ris1,pos_ue(:,userIdx),v,v) + ynlos(:,i,userIdx)+chanRISToUE(x_ris_out2(:,i,userIdx),pos_ris2,pos_ue(:,userIdx),v,v);
        SNRlos2ris(i,userIdx) = pow2db(bandpower(ylos2ris(:,i,userIdx))) - N0dB;
    end
end

% Calculate the average SNR for each user with 1 ris
averageSNRUsers1 = mean(SNRlos1ris, 1);

% Display the average SNR for each user
for userIdx = 1:numUsers
    disp(['Average SNR for User ', num2str(userIdx), ' with 1 RIS: ', num2str(averageSNRUsers1(userIdx)), ' dB']);
end

% Calculate the average SNR for each user with 2 ris
averageSNRUsers2 = mean(SNRlos2ris, 1);

% Display the average SNR for each user
for userIdx = 1:numUsers
    disp(['Average SNR for User ', num2str(userIdx), ' with 2 RIS: ', num2str(averageSNRUsers2(userIdx)), ' dB']);
end


%% Plot SNR with RIS vs. User Locations

% Extract user positions and SNR values with RIS
user_positions = pos_ue'; % Transpose to match user dimension
snr_with_1ris = averageSNRUsers1; % Average SNR with 1 RIS for each user
snr_with_2ris= averageSNRUsers2;% Average SNR with 2 RIS for each user
% Plot
figure;
scatter(user_positions(:,1), snr_with_1ris, 'filled');
xlabel('Distance (meters)');
ylabel('SNR with RIS (dB)');
title('SNR with RIS vs. User Location');
grid on;
hold on
scatter(user_positions(:,1), snr_with_2ris, 'filled');
%% Plot Average SNR vs. Distance and SNR with RIS vs. User Locations

% Extract user positions and SNR values with RIS
user_positions = pos_ue'; % Transpose to match user dimension
snr_with_1ris = averageSNRUsers1; % Average SNR with RIS for each user

% Create a figure with 2 subplots
figure;

% Plot 1: Average SNR vs. Distance
subplot(1, 3, 1);
scatter(distances, average_snr_dB, 'filled');
xlabel('Distance (meters)');
ylabel('Average SNR (dB)');
title('Average SNR vs. Distance');
grid on;

% Plot 2: SNR with 1 RIS vs. User Locations
subplot(1, 3, 2);
scatter(user_positions(:,1), snr_with_1ris, 'filled');
xlabel('Distance (meters)');
ylabel('SNR with one RIS (dB)');
title('SNR with RIS vs. User Location');
grid on;



% Plot 3: SNR with 1 RIS vs. User Locations
subplot(1, 3, 3);
scatter(user_positions(:,1), snr_with_2ris, 'filled');
xlabel('Distance (meters)');
ylabel('SNR with two RIS (dB)');
title('SNR with RIS vs. User Location');
grid on;

% Adjust layout to make room for labels
sgtitle('SNR Analysis');
