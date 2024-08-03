clear all
close all
clc

%% Definition of system parameters
numTx = 2;        % Number of transmit antennas
numRx = 2;        % Number of receive antennas
numUsers = 100;     % Number of users
fc = 5e9;         % Carrier frequency
lambda = physconst('LightSpeed') / fc;
distances = 100:2:299;  % Generate random distances between 100 and 300 meters
numDistances=length(distances);
meanShadowing_dB = 0;      % Mean shadowing in dB
stdShadowing_dB = 4;       % Standard deviation of shadowing in dB
numRealizations = 5000;    % Number of Monte Carlo realizations

%% Desired Transmit Power
desiredPower_dBm = 46; % Desired transmit power in dBm
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
        
        % Transmited signal
        xt = xt1*sqrt(desiredPower_Watts);
        numbSymbols = length(xt(:,1));
        bandpower(xt)

        %% Path Loss Initialization
        pathLoss_dB = fspl(distances, lambda);
        tx_pathloss = zeros(numbSymbols, numTx, numDistances);
        for disIdx = 1:numDistances
            pathLossFactor = 10^(-pathLoss_dB(disIdx) / 10);
            for tx = 1:numTx
                tx_pathloss(:, tx, disIdx) = xt(:, tx) * sqrt(pathLossFactor);
            end
        end

        %% Shadowing
        shadowing_dB = stdShadowing_dB * randn(numRealizations, numDistances, numTx) + meanShadowing_dB;
        tx_shadowing = zeros(numbSymbols, numTx, numDistances);
        for disIdx = 1:numDistances
            for tx = 1:numTx
                for realization = 1:numRealizations
                    shadowingFactor = 10^(-shadowing_dB(realization, disIdx, tx) / 20);
                    shadowing_effect = tx_pathloss(:, tx, disIdx) * sqrt(shadowingFactor);
                    tx_shadowing(:, tx, disIdx) = tx_shadowing(:, tx, disIdx) + shadowing_effect;
                end
                tx_shadowing(:, tx, disIdx) = tx_shadowing(:, tx, disIdx) / numRealizations;
            end
        end

        %% Apply Rayleigh Fading
        received_signal = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            faded_signal = zeros(numbSymbols, numRx);
             for realization = 1:numRealizations
                 % Generate Rayleigh fading matrix for this realization
                 rayleigh_fading = (randn(numRx, numTx) + 1j * randn(numRx, numTx)) / sqrt(2);
                 % Apply Rayleigh fading to the signal
                   for rx = 1:numRx
                      for tx = 1:numTx
                        faded_signal(:, rx) = faded_signal(:, rx) + rayleigh_fading(rx, tx) * tx_shadowing(:, tx, disIdx);
                      end
                   end
             end    
            faded_signal = faded_signal / numRealizations;
            received_signal(:, :, disIdx) = faded_signal;
        end

        %% Add AWGN
        noise_variance = db2pow(-120);
        awgn_noise = sqrt(noise_variance / 2) * (randn(numbSymbols, numRx, numDistances) + 1j * randn(numbSymbols, numRx, numDistances));

        %% Calculate SNR
        snr_dB = zeros(numRx, numDistances);
        received_with_noise = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            for rx = 1:numRx
                received_with_noise(:, rx, disIdx) = received_signal(:, rx, disIdx) + awgn_noise(:, rx, disIdx);
                signal_power_Watts = bandpower(received_signal(:, rx, disIdx));
                noise_power_Watts = noise_variance;
                snr_dB(rx, disIdx) = 10 * log10(signal_power_Watts / noise_power_Watts);
            end
        end
        average_snr_dB = mean(snr_dB, 1);

        %% Demodulate and Calculate BER
        demodulated_bits = demodulation(received_with_noise, modScheme, M);
        
        ber = zeros(numDistances, 1);
        throughput = zeros(numDistances, 1);
        tx_bits = bit_stream_tx(:, :);

        for disIdx = 1:numDistances
            demod_bits = demodulated_bits(:, :, disIdx);
            bit_errors = sum(tx_bits ~= demod_bits, 'all');
            total_bits = numel(tx_bits);
            ber(disIdx) = bit_errors / total_bits;
            data_rate = log2(M);
            throughput(disIdx) = data_rate * (1 - ber(disIdx));
        end

        %% Store Results
        results.(modScheme).(['M' num2str(M)]) = struct('SNR', average_snr_dB, 'BER', ber, 'Throughput', throughput);
        %% Plot SNR vs Distance
        snr_values = zeros(numDistances, 1);
        for disIdx = 1:numDistances
            snr_values(disIdx) = mean(results.(modScheme).(['M' num2str(M)]).SNR(:, disIdx));
        end
        % Interpolation to get smooth curve
        interpDistances = linspace(min(distances), max(distances), 100);  % 100 points for smooth curve
        interpoMethod='pchip';
        interpSNR = interp1(distances, snr_values, interpDistances, interpMethod);
        subplot(length(modulationSchemes), length(orders), (schemeIdx-1)*length(orders) + orderIdx);
        plot(distances, snr_values, 'o', interpDistances, interpSNR, '-'); 
        title([modScheme ' M' num2str(M)]);
        xlabel('Distance (m)');
        ylabel('SNR (dB)');
        legend('Original', 'Interpolated');
        grid on;


        
    end
end

sgtitle('SNR vs Distance for Different Modulation Schemes and Orders');
