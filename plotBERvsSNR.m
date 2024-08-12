function plotBERvsSNR(results, distances)
    figure;
    title('BER vs. SNR for Different Modulation Schemes and Orders');
    xlabel('SNR (dB)');
    ylabel('BER');

    % Loop through modulation schemes, orders, and coding rates
    modulationSchemes = fieldnames(results);
    colors = lines(length(modulationSchemes)); % Different colors for each scheme

    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        orders = fieldnames(results.(modScheme));
        
        for orderIdx = 1:length(orders)
            orderFieldName = orders{orderIdx};
            codingRates = fieldnames(results.(modScheme).(orderFieldName));

            for rateIdx = 1:length(codingRates)
                codingRateFieldName = codingRates{rateIdx};
                ber = results.(modScheme).(orderFieldName).(codingRateFieldName).BER_1RIS;
                snr_dB = results.(modScheme).(orderFieldName).(codingRateFieldName).SNR_1RIS;

                % Plot BER vs SNR
                plot(snr_dB, ber, 'o-', 'DisplayName', sprintf('%s M%d Rate %.0f%%', modScheme, str2double(orderFieldName(2:end)), str2double(codingRateFieldName(5:end))));
            end
        end
    end

    legend('show');
    grid on;
    hold off;
end
