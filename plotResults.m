function plotResults(results, distances,codingRates)
% Extract modulation schemes
    modulationSchemes = fieldnames(results);
    interpMethod = 'pchip';  % Define interpolation method
    interpDistances = linspace(min(distances), max(distances), 100);  % 100 points for smooth curves
    
    % Plot for different modulation orders with same coding rate
    figure;
    numModSchemes = length(modulationSchemes);
    numCodingRates = length(codingRates);

    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        modulationOrders = fieldnames(results.(modScheme));
        
        % Iterate through each coding rate
        for rateIdx = 1:numCodingRates
            codingRate = codingRates(rateIdx);
            codingRateFieldName = sprintf('codingRate%.0f', codingRate * 100);
            subplot(numModSchemes, numCodingRates, (schemeIdx-1) * numCodingRates + rateIdx);
            hold on;
            % Iterate through each modulation order
            for orderIdx = 1:length(modulationOrders)
                M = modulationOrders{orderIdx};
                snr_noRIS = results.(modScheme).(M).(codingRateFieldName).SNR_NoRIS;
                snr_1RIS = results.(modScheme).(M).(codingRateFieldName).SNR_1RIS;
                snr_2RIS = results.(modScheme).(M).(codingRateFieldName).SNR_2RIS;
                
                % Interpolation
                interpSNR_noRIS = interp1(distances, snr_noRIS, interpDistances, interpMethod);
                interpSNR_1RIS = interp1(distances, snr_1RIS, interpDistances, interpMethod);
                interpSNR_2RIS = interp1(distances, snr_2RIS, interpDistances, interpMethod);
                
                % Plot
                plot(interpDistances, interpSNR_noRIS, '-');
                plot(interpDistances, interpSNR_1RIS, '--');
                plot(interpDistances, interpSNR_2RIS, ':');
            end
            
            title(['Coding Rate ', num2str(codingRates(rateIdx))]);
            xlabel('Distance (m)');
            ylabel('SNR (dB)');
            legend(modulationOrders, 'Location', 'best');
            grid on;
        end
    end
    
    sgtitle('SNR vs Distance for Different Modulation Orders and Coding Rates');
    
    % Additional plots can be created similarly, e.g., for different coding rates with the same modulation order
end