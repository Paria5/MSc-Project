clear all
close all
clc
%% Definition of system parameters
numTx = 2;        % Number of transmit antennas
numRx = 2;        % Number of receive antennas
numUsers = 100;     % Number of users
fc = 5e9;         % Carrier frequency
lambda = physconst('LightSpeed') / fc;
distances = randi([100, 300], numUsers, 1);  % Generate random distances between 100 and 300 meters
meanShadowing_dB = 0;      % Mean shadowing in dB
stdShadowing_dB = 0;       % Standard deviation of shadowing in dB
numRealizations = 1000;    % Number of Monte Carlo realizations
%% Desired Transmit Power
desiredPower_dBm = 46; % Desired transmit power in dBm
desiredPower_Watts = 10^((desiredPower_dBm - 30) / 10); % Convert dBm to Watts
%% Bit Generation
numBits=1000;
bit_stream_tx = randi([0 1], numTx, numBits);
M=2; %bpsk
k=log2(M);
numbSymbols=numBits/k;
modulated_tx =pskmod(bit_stream_tx,M)'; 
% signal
tx = phased.Transmitter('PeakPower',1,'Gain',3);
xt = tx(modulated_tx);
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
%% Apply Rayleigh Fading
received_signal = zeros(numbSymbols,numRx, numUsers);
numRealization_r=10000;
for userIdx = 1:numUsers
    % Apply Rayleigh fading
    faded_signal = zeros(numbSymbols,numRx);
    for realization = 1:numRealization_r
        % Generate Rayleigh fading matrix for this realization
        rayleigh_fading = (randn(numRx, numTx) + 1j * randn(numRx, numTx)) / sqrt(2);
        
        % Apply Rayleigh fading to the signal
        for rx = 1:numRx
            for tx = 1:numTx
                faded_signal(:, rx) = faded_signal(:, rx) + rayleigh_fading(rx, tx) * tx_shadowing(:, tx, userIdx);
            end
        end
    end    
    % Average the faded signal over all realizations
    faded_signal = faded_signal / numRealization_r;
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
noise_variance = db2pow(-100);
awgn_noise = sqrt(noise_variance / 2) * (randn(numbSymbols, numRx, numUsers) + 1j * randn(numbSymbols, numRx, numUsers));

%% Calculate Noise Power in dBm
noise_power_Watts = db2pow(-100);
%noise_power_dBm = -100;
%% Calculate SNR
snr_dB = zeros(numRx, numUsers);
received_with_noise=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for rx = 1:numRx
        % Calculate noise power
        %noise_power_Watts = noise_variance;
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

dv = -10;
du = sqrt((distances.^2)-dv^2);
pos_ue = [du.'; dv*ones(1,numUsers); zeros(1,numUsers)];

% scene
dbr = 100;
pos_ap = [0;0;0];
pos_ris = [dbr;0;0]; 
v = zeros(3,1);

% compute the range and angle of the RIS from the base station and the UE
[r_ap_ris,ang_ap_ris] = rangeangle(pos_ap,pos_ris);
[r_ue_ris,ang_ue_ris] = rangeangle(pos_ue,pos_ris);

fs = 10e6;
c = physconst('lightspeed');

% channel
chanAPToRIS = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);
chanRISToUE = phased.FreeSpace('SampleRate',fs,'PropagationSpeed',c,'MaximumDistanceSource','Property','MaximumDistance',500);

N0dB=-4;
rcoeff_ris = ones(Nr*Nc,1);
x_ris_in=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_in(:,i,userIdx) = chanAPToRIS(xt(:,i),pos_ap,pos_ris,v,v);
    end
end
x_ris_out=zeros(numbSymbols,numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        x_ris_out(:,i,userIdx) = ris(x_ris_in(:,i,userIdx),ang_ap_ris,ang_ue_ris(:,userIdx),rcoeff_ris);
    end
end
ylosris=zeros(numbSymbols,numRx,numUsers);
ynlos=received_signal;
SNRlosris=zeros(numRx,numUsers);
for userIdx = 1:numUsers
    for i=1:numRx
        ylosris(:,i,userIdx) = chanRISToUE(x_ris_out(:,i,userIdx),pos_ris,pos_ue(:,userIdx),v,v) + ynlos(:,i,userIdx);
        SNRlosris(i,userIdx) = pow2db(bandpower(ylosris(:,i,userIdx))) - N0dB;
    end
end
% Calculate the average SNR for each user
averageSNRUsers = mean(SNRlosris, 1);

% Display the average SNR for each user
for userIdx = 1:numUsers
    disp(['Average SNR for User ', num2str(userIdx), ' with RIS: ', num2str(averageSNRUsers(userIdx)), ' dB']);
end
%% Plot SNR with RIS vs. User Locations

% Extract user positions and SNR values with RIS
user_positions = pos_ue'; % Transpose to match user dimension
snr_with_ris = averageSNRUsers; % Average SNR with RIS for each user

% Plot
figure;
scatter(user_positions(:,1), snr_with_ris, 'filled');
xlabel('Distance (meters)');
ylabel('SNR with RIS (dB)');
title('SNR with RIS vs. User Location');
grid on;
%% Plot Average SNR vs. Distance and SNR with RIS vs. User Locations

% Extract user positions and SNR values with RIS
user_positions = pos_ue'; % Transpose to match user dimension
snr_with_ris = averageSNRUsers; % Average SNR with RIS for each user

% Create a figure with 2 subplots
figure;

% Plot 1: Average SNR vs. Distance
subplot(1, 2, 1);
scatter(distances, average_snr_dB, 'filled');
xlabel('Distance (meters)');
ylabel('Average SNR (dB)');
title('Average SNR vs. Distance');
grid on;

% Plot 2: SNR with RIS vs. User Locations
subplot(1, 2, 2);
scatter(user_positions(:,1), snr_with_ris, 'filled');
xlabel('Distance (meters)');
ylabel('SNR with RIS (dB)');
title('SNR with RIS vs. User Location');
grid on;

% Adjust layout to make room for labels
sgtitle('SNR Analysis');
