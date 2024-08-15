function generateSNRTable(results, allSNRsorted, outputFileName)
    % Generate MCS index mapping
    index = 1;
    numMCS =15;
    rows=length(allSNRsorted)/numMCS;
    modulationSchemes = fieldnames(results);  % Get modulation schemes

    mcsIndices = struct(); % to store the index mapping

    for schemeIdx = 1:length(modulationSchemes)
        modScheme = modulationSchemes{schemeIdx};
        modulationOrders = fieldnames(results.(modScheme));  % Get modulation orders
        
        for orderIdx = 1:length(modulationOrders)
            M = modulationOrders{orderIdx};
            codingRates = fieldnames(results.(modScheme).(M));  % Get coding rates
            for rateIdx = 1:length(codingRates)
                codingRate = codingRates(rateIdx);
                mcsIndices.(sprintf('index%d', index)) = struct('Modulation', modScheme, 'Order', M, 'CodingRate', codingRate);
                index = index + 1;
            end
        end
    end

    % Initialize SNR results storage
    SNR_table = table('Size', [length(allSNRsorted)/numMCS, numMCS], ...
        'VariableTypes', repmat({'double'}, 1, numMCS), ...
        'VariableNames', arrayfun(@(i) num2str(i), 1:numMCS, 'UniformOutput', false))
    
% Fill the table with SNR values
    for mcsIdx = 1:numMCS
        % Determine the indices for this MCS column
        startIdx = (mcsIdx - 1) * (rows) + 1;
        endIdx = mcsIdx * (rows);
        
        % Extract values for this MCS index
        snrValues = allSNRsorted(startIdx:endIdx);
        
        % Assign to the corresponding column in the table
        SNR_table{:, mcsIdx} = snrValues.';
    end

    % Write table to file
    writetable(SNR_table, outputFileName);
end
