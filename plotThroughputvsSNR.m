function plotThroughputvsSNR(results, distances)
    figure;
    hold on;
    title('Throughput vs. SNR for Different Modulation Schemes and Orders');
    xlabel('SNR (dB)');
    ylabel('Throughput (bps/Hz)');

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
                throughput = results.(modScheme).(orderFieldName).(codingRateFieldName).Throughput_NoRIS;
                snr_dB = results.(modScheme).(orderFieldName).(codingRateFieldName).SNR_NoRIS;

                % Plot Throughput vs SNR
                plot(snr_dB, throughput, 'o-', 'DisplayName', sprintf('%s M%d Rate %.0f%%', modScheme, str2double(orderFieldName(2:end)), str2double(codingRateFieldName(5:end))));
            end
        end
    end

    legend('show');
    grid on;
    hold off;
end
