function generateSNRTable(results, distances, modulationSchemes, modulationOrders, codingRates, outputFileName)
    % Generate MCS index mapping
    index = 1;
    numMCS = 24;
    mcsIndices = struct(); % to store the index mapping

    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        orders = modulationOrders.(modScheme);
        
        for orderIdx = 1:length(orders)
            M = orders(orderIdx);
            for rateIdx = 1:length(codingRates)
                codingRate = codingRates(rateIdx);
                mcsIndices.(sprintf('Index%d', index)) = struct('Modulation', modScheme, 'Order', M, 'CodingRate', codingRate);
                index = index + 1;
            end
        end
    end

    % Initialize SNR results storage
    numDistances = length(distances);
    SNR_table = table('Size', [numDistances, numMCS], ...
        'VariableTypes', repmat({'double'}, 1, numMCS), ...
        'VariableNames', arrayfun(@(i) sprintf('Index%d', i), 1:numMCS, 'UniformOutput', false));
    
    % Populate the table with SNR values
    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        orders = modulationOrders.(modScheme);
        
        for orderIdx = 1:length(orders)
            M = orders(orderIdx);
            
            for rateIdx = 1:length(codingRates)
                codingRate = codingRates(rateIdx);
                
                % Find the index for the current MCS
                indexName = sprintf('Index%d', find(arrayfun(@(i) isequal(mcsIndices.(sprintf('Index%d', i)).Modulation, modScheme) && ...
                                                                 mcsIndices.(sprintf('Index%d', i)).Order == M && ...
                                                                 mcsIndices.(sprintf('Index%d', i)).CodingRate == codingRate, ...
                                                1:numMCS)));
                
                % Get the results for this MCS
                codingRateFieldName = sprintf('codingRate%.0f', codingRate * 100);
                
                if isfield(results.(modScheme).(['M' num2str(M)]), codingRateFieldName)
                    snrData = results.(modScheme).(['M' num2str(M)]).(codingRateFieldName).SNR_NoRIS; % Adjust if RIS is used
                    % Store SNR values in the table
                    SNR_table.(indexName) = num2cell(snrData');
                else
                    warning('Field %s does not exist in results for Modulation %s, Order %d, CodingRate %.2f', ...
                        codingRateFieldName, modScheme, M, codingRate);
                end
            end
        end
    end
    
    % Save the results to a file
    writetable(SNR_table, outputFileName);
end
