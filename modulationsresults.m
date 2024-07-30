clear all
close all
clc

%% Definition of system parameters
numTx = 2;        % Number of transmit antennas
numRx = 2;        % Number of receive antennas
numUsers = 50;     % Number of users
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
numBits = 1000;
bit_stream_tx = randi([0 1], numBits, numTx);

% Define modulation schemes and their corresponding orders
modulationSchemes = {'qam', 'psk'};
modulationOrders = struct('qam', [4, 16, 256], 'psk', [2, 4, 16]);

% Initialize results storage
results = struct();

for schemeIdx = 1:length(modulationSchemes)
    modScheme = modulationSchemes{schemeIdx};
    orders = modulationOrders.(modScheme);
    
    for orderIdx = 1:length(orders)
        M = orders(orderIdx);
        
        % Reshape the bit stream to fit the modulator
        bit_stream_tx_reshaped = reshape(bit_stream_tx, [], numTx);
        
        % Modulate the bit stream
        modulated_tx = [];
        for tx = 1:numTx
            modulated_tx = [modulated_tx, modulation(bit_stream_tx_reshaped(:, tx), modScheme, M)];
        end
        
        % Normalize the modulated signal
        modulated_power = mean(abs(modulated_tx).^2);
        xt1 = modulated_tx ./ sqrt(modulated_power);
        
        % Transmit signal
        tx = phased.Transmitter('PeakPower',1,'Gain',3);
        xt = tx(xt1);
        numbSymbols = length(xt(:,1));

        %% Path Loss Initialization
        pathLoss_dB = fspl(distances, lambda);
        tx_pathloss = zeros(numbSymbols, numTx, numUsers);
        for userIdx = 1:numUsers
            pathLossFactor = 10^(-pathLoss_dB(userIdx) / 10);
            for tx = 1:numTx
                tx_pathloss(:, tx, userIdx) = xt(:, tx) * sqrt(pathLossFactor);
            end
        end

        %% Shadowing
        shadowing_dB = stdShadowing_dB * randn(numRealizations, numUsers, numTx) + meanShadowing_dB;
        tx_shadowing = zeros(numbSymbols, numTx, numUsers);
        for userIdx = 1:numUsers
            for tx = 1:numTx
                for realization = 1:numRealizations
                    shadowingFactor = 10^(-shadowing_dB(realization, userIdx, tx) / 20);
                    shadowing_effect = tx_pathloss(:, tx, userIdx) * sqrt(shadowingFactor);
                    tx_shadowing(:, tx, userIdx) = tx_shadowing(:, tx, userIdx) + shadowing_effect;
                end
                tx_shadowing(:, tx, userIdx) = tx_shadowing(:, tx, userIdx) / numRealizations;
            end
        end

        %% Rayleigh Fading
        rayleigh_fading = (randn(numRx, numTx) + 1j * randn(numRx, numTx)) / sqrt(2);

        %% Apply Rayleigh Fading
        received_signal = zeros(numbSymbols, numRx, numUsers);
        for userIdx = 1:numUsers
            faded_signal = zeros(numbSymbols, numRx);
            for rx = 1:numRx
                for tx = 1:numTx
                    faded_signal(:, rx) = faded_signal(:, rx) + rayleigh_fading(rx, tx) * tx_shadowing(:, tx, userIdx);
                end
            end
            received_signal(:, :, userIdx) = faded_signal;
        end

        %% Add AWGN
        noise_variance = 0.5;
        awgn_noise = sqrt(noise_variance / 2) * (randn(numbSymbols, numRx, numUsers) + 1j * randn(numbSymbols, numRx, numUsers));

        %% Calculate SNR
        snr_dB = zeros(numRx, numUsers);
        received_with_noise = zeros(numbSymbols, numRx, numUsers);
        for userIdx = 1:numUsers
            for rx = 1:numRx
                received_with_noise(:, rx, userIdx) = received_signal(:, rx, userIdx) + awgn_noise(:, rx, userIdx);
                signal_power_Watts = bandpower(received_signal(:, rx, userIdx));
                noise_power_Watts = noise_variance;
                snr_dB(rx, userIdx) = 10 * log10(signal_power_Watts / noise_power_Watts);
            end
        end
        average_snr_dB = mean(snr_dB, 1);

        %% Demodulate and Calculate BER
        demodulated_bits = demodulation(received_with_noise, modScheme, M);
        
        ber = zeros(numUsers, 1);
        throughput = zeros(numUsers, 1);
        tx_bits = bit_stream_tx(:, :);

        for userIdx = 1:numUsers
            demod_bits = demodulated_bits(:, :, userIdx);
            bit_errors = sum(tx_bits ~= demod_bits, 'all');
            total_bits = numel(tx_bits);
            ber(userIdx) = bit_errors / total_bits;
            data_rate = log2(M);
            throughput(userIdx) = data_rate * (1 - ber(userIdx));
        end

        %% Store Results
        results.(modScheme).(['M' num2str(M)]) = struct('SNR', average_snr_dB, 'BER', ber, 'Throughput', throughput);
    end
end
%% Plotting SNR vs Distance for Both QAM and PSK on the Same Graph
figure;
hold on;
markers = {'o', 's', 'd', 'x', '^', 'v', '<', '>', 'p', 'h'};
lineStyles = {'-', '--', '-.', ':'};  % Define line styles
markerIdx = 1;
lineStyleIdx = 1;

for schemeIdx = 1:length(modulationSchemes)
    modScheme = modulationSchemes{schemeIdx};
    orders = modulationOrders.(modScheme);
    
    for orderIdx = 1:length(orders)
        M = orders(orderIdx);
        
        % Extract SNR values for the current modulation scheme and order
        snr_values = zeros(numUsers, 1);
        for userIdx = 1:numUsers
            snr_values(userIdx) = mean(results.(modScheme).(['M' num2str(M)]).SNR(:, userIdx));
        end
        
        % Sort distances and corresponding average SNR values
        [uniqueDistances, sortIdx] = unique(distances);  % Ensure unique distances
        avgSNR = zeros(length(uniqueDistances), 1);
        
        for i = 1:length(uniqueDistances)
            % Find all indices corresponding to this distance
            indices = (distances == uniqueDistances(i));
            avgSNR(i) = mean(snr_values(indices));
        end
        
        % Interpolation to get smooth curve
        interpDistances = linspace(min(uniqueDistances), max(uniqueDistances), 100);  % 100 points for smooth curve
        interpSNR = interp1(uniqueDistances, avgSNR, interpDistances, 'pchip', 'extrap');  % 'pchip' for piecewise cubic Hermite interpolation
        
        % Plot with lines connecting the points
        plot(interpDistances, interpSNR, ...
             'DisplayName', [upper(modScheme) ', M=' num2str(M)], ...
             'LineStyle', lineStyles{lineStyleIdx}, ...
             'LineWidth', 1.5, ...
             'Marker', markers{markerIdx}, ...
             'MarkerSize', 8, ...
             'MarkerIndices', 1:length(uniqueDistances));  % Only mark the original points

        markerIdx = markerIdx + 1;
        if markerIdx > length(markers)
            markerIdx = 1;
            lineStyleIdx = lineStyleIdx + 1;
            if lineStyleIdx > length(lineStyles)
                lineStyleIdx = 1;  % Reset line styles if needed
            end
        end
    end
end

xlabel('Distance (m)');
ylabel('Average SNR (dB)');
title('Average SNR vs Distance for QAM and PSK');
legend show;
grid on;
hold off;