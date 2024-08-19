function allSNRsSorted = extractAndSortSNRwithRIS(results)
    allSNRs = [];  % Initialize an empty array to store SNR values

    modulationSchemes = fieldnames(results);  % Get modulation schemes
    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        
        modulationOrders = fieldnames(results.(modScheme));  % Get modulation orders
        for orderIdx = 1:length(modulationOrders)
            modOrder = modulationOrders{orderIdx};
            
            codingRates = fieldnames(results.(modScheme).(modOrder));  % Get coding rates
            for rateIdx = 1:length(codingRates)
                codingRate = codingRates{rateIdx};
                
                % Extract SNR values without RIS
                snrValues = results.(modScheme).(modOrder).(codingRate).SNR_1RIS;
                
                % Append these SNR values to the allSNRs array
                allSNRs = [allSNRs, snrValues];
            end
        end
    end
    
    % Sort all SNRs in descending order
    allSNRsSorted = sort(allSNRs, 'ascend');
end