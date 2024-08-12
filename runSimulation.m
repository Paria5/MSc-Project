function results = runSimulation(numUsers,numTx,ris1x,ris1y, ris1z, ris2x, ris2y,ris2z)
    % Main simulation code
    clearvars -except numTx numUsers ris1x ris1y ris1z ris2x ris2y ris2z
    close all
    clc
    fc = 5e9;         % Carrier frequency
    numRx=numTx;
    pos_ris1 = [ris1x;ris1y;ris1z]; 
    pos_ris2 = [ris2x;ris2y;ris2z];
    fs = 10e6;
c = physconst('lightspeed');
lambda = c/ fc;
distances = 100:20:299;  % Generate random distances between 100 and 300 meters
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
modulationOrders = struct('qam', [4, 16,64, 256], 'psk', [2, 4,8, 16]);
codingRates = [1/2, 2/3, 3/4];

% Initialize results storage
results = struct();

for schemeIdx = 1:length(modulationSchemes)
    modScheme = modulationSchemes{schemeIdx};
    orders = modulationOrders.(modScheme);
    
    for orderIdx = 1:length(orders)
        M = orders(orderIdx);
        
        %for each coding rate
     for rateIdx=1:length(codingRates)
           codingRate=codingRates(rateIdx);
            
        % Reshape the bit stream to fit the modulator
        bit_stream_tx_reshaped = reshape(bit_stream_tx, [], numTx);
        
        % Modulate the bit stream
        modulated_tx = [];
        for tx = 1:numTx
            modulated_tx = [modulated_tx, modulation(bit_stream_tx_reshaped(:, tx), modScheme, M,codingRate)];
        end
        
        % Normalize the modulated signal
        modulated_power = mean(abs(modulated_tx).^2);
        xt1 = modulated_tx ./ sqrt(modulated_power);
        
        % Transmited signal
        xt = xt1*sqrt(desiredPower_Watts/numUsers);
        transmitted_pow_PU=bandpower(xt);
        numbSymbols = length(xt(:,1));

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
        demodulated_bits = demodulation(received_with_noise, modScheme, M,numBits,codingRate);
       
        ber = zeros(numDistances, 1);
        throughput = zeros(numDistances, 1);
        tx_bits = bit_stream_tx(:, :);

        for disIdx = 1:numDistances
            demod_bits = demodulated_bits(:, :, disIdx);
            bit_errors = sum(tx_bits ~= demod_bits, 'all');
            total_bits = numel(tx_bits);
            ber(disIdx) = bit_errors / total_bits;
            data_rate = log2(M);
            % Effective data rate considering the coding rate
            effective_data_rate = data_rate * (1 - codingRate);
            throughput(disIdx) = effective_data_rate * (1 - ber(disIdx));
        end

        %% Store Results without RIS
        codingRateFieldName = sprintf('codingRate%.0f', codingRate * 100); % Format codingRate as a valid field name
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).SNR_NoRIS = average_snr_dB;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).BER = ber;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).Throughput = throughput;

        %% RIS Setup
        Nr = 10;
        Nc = 20;
        dr = 0.5 * lambda;
        dc = 0.5 * lambda;

        % Construct surface
        ris = helperRISSurface('Size', [Nr Nc], 'ElementSpacing', [dr dc], ...
            'ReflectorElement', phased.IsotropicAntennaElement, 'OperatingFrequency', fc);

        dv = 0;
        du = sqrt((distances.^2) - dv^2);
        pos_ue = [du; dv*ones(1, numDistances); zeros(1, numDistances)];

        % Scene for RIS1
        %dbr1 = 100; 
        pos_ap = [0; 0; 0];
        v = zeros(3, 1);

        % Compute the range and angle of the RIS from the base station and the UE
        [r_ap_ris, ang_ap_ris] = rangeangle(pos_ap, pos_ris1);
        [r_ue_ris, ang_ue_ris] = rangeangle(pos_ue, pos_ris1);

        % Channel
        chanAPToRIS = phased.FreeSpace('SampleRate', fs, 'PropagationSpeed', c, 'MaximumDistanceSource', 'Property', 'MaximumDistance', 500);
        chanRISToUE = phased.FreeSpace('SampleRate', fs, 'PropagationSpeed', c, 'MaximumDistanceSource', 'Property', 'MaximumDistance', 500);

        N0dB = -120;
        rcoeff_ris = ones(Nr * Nc, 1);
        x_ris_in1 = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                x_ris_in1(:, i, disIdx) = chanAPToRIS(xt(:, i), pos_ap, pos_ris1, v, v);
            end
        end
        x_ris_out1 = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                x_ris_out1(:, i, disIdx) = ris(x_ris_in1(:, i, disIdx), ang_ap_ris, ang_ue_ris(:, disIdx), rcoeff_ris);
            end
        end
        ynlos = received_signal;

        ylos1ris = zeros(numbSymbols, numRx, numDistances);
        SNRlos1ris = zeros(numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                ylos1ris(:, i, disIdx) = chanRISToUE(x_ris_out1(:, i, disIdx), pos_ris1, pos_ue(:, disIdx), v, v) + ynlos(:, i, disIdx);
                SNRlos1ris(i, disIdx) = pow2db(bandpower(ylos1ris(:, i, disIdx))) - N0dB;
            end
        end
        meanSNR_with1RIS = mean(SNRlos1ris, 1);
        
        %% Calculate BER and Throughput for 1 RIS
        demodulated_bits_1RIS = demodulation(ylos1ris, modScheme, M,numBits,codingRate);
        
        ber_1RIS = zeros(numDistances, 1);
        throughput_1RIS = zeros(numDistances, 1);

        for disIdx = 1:numDistances
            demod_bits_1RIS = demodulated_bits_1RIS(:, :, disIdx);
            bit_errors_1RIS = sum(tx_bits ~= demod_bits_1RIS, 'all');
            total_bits = numel(tx_bits);
            ber_1RIS(disIdx) = bit_errors_1RIS / total_bits;
            throughput_1RIS(disIdx) = data_rate * (1 - ber_1RIS(disIdx));
        end

        %% Store Results with 1 RIS
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).SNR_1RIS = meanSNR_with1RIS;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).BER_1RIS = ber_1RIS;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).Throughput_1RIS = throughput_1RIS;

        %% Add one more RIS
        x_ris_in2 = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                x_ris_in2(:, i, disIdx) = chanAPToRIS(xt(:, i), pos_ap, pos_ris2, v, v);
            end
        end
        x_ris_out2 = zeros(numbSymbols, numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                x_ris_out2(:, i, disIdx) = ris(x_ris_in2(:, i, disIdx), ang_ap_ris, ang_ue_ris(:, disIdx), rcoeff_ris);
            end
        end
        ylos2ris = zeros(numbSymbols, numRx, numDistances);
        SNRlos2ris = zeros(numRx, numDistances);
        for disIdx = 1:numDistances
            for i = 1:numRx
                ylos2ris(:, i, disIdx) = chanRISToUE(x_ris_out1(:, i, disIdx), pos_ris1, pos_ue(:, disIdx), v, v) + ynlos(:, i, disIdx) + chanRISToUE(x_ris_out2(:, i, disIdx), pos_ris2, pos_ue(:, disIdx), v, v);
                SNRlos2ris(i, disIdx) = pow2db(bandpower(ylos2ris(:, i, disIdx))) - N0dB;
            end
        end
        meanSNR_with2RIS = mean(SNRlos2ris, 1);
        %% Calculate BER and Throughput for 2 RIS
        demodulated_bits_2RIS = demodulation(ylos2ris, modScheme, M,numBits,codingRate);
        
        ber_2RIS = zeros(numDistances, 1);
        throughput_2RIS = zeros(numDistances, 1);

        for disIdx = 1:numDistances
            demod_bits_2RIS = demodulated_bits_2RIS(:, :, disIdx);
            bit_errors_2RIS = sum(tx_bits ~= demod_bits_2RIS, 'all');
            total_bits = numel(tx_bits);
            ber_2RIS(disIdx) = bit_errors_2RIS / total_bits;
            throughput_2RIS(disIdx) = data_rate * (1 - ber_2RIS(disIdx));
        end

        %% Store Results with 2 RIS
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).SNR_2RIS = meanSNR_with2RIS;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).BER_2RIS = ber_2RIS;
        results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).Throughput_2RIS = throughput_2RIS;
    end
    end
end
save('simulation_results.mat', 'results');
plotResults(results,distances,codingRates);
sgtitle('SNR vs Distance for Different Modulation Schemes and Orders');
end