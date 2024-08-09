function bit_stream = demodulation(modulated_signal, modulation_type, M,numBitsOriginal)
    [numSymbols, numTx, numDistance] = size(modulated_signal);
    bitsPerSymbol = log2(M);
    
    % Validate inputs
    validateattributes(modulation_type, {'char', 'string'}, {'nonempty'});
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'positive'});
    
    % Initialize matrix for demodulated bits
    numBits = numSymbols * bitsPerSymbol;
    bit_stream = zeros(numBitsOriginal, numTx, numDistance); 

    % Process each user's signal
    for disIdx = 1:numDistance
        for tx = 1:numTx
            % Extract the current column (transmitter) for the current user
            current_modulated_signal = modulated_signal(:, tx, disIdx);

            % Perform demodulation based on type
            if strcmp(modulation_type, 'psk') % PSK demodulation
                dataStream = pskdemod(current_modulated_signal, M);
            elseif strcmp(modulation_type, 'qam') % QAM demodulation
                dataStream = qamdemod(current_modulated_signal, M);
            else
                error('Unsupported modulation type.');
            end

            % Convert symbols to bits
            bitStreamArray = de2bi(dataStream, bitsPerSymbol, 'left-msb').';
            bitStreamArray = bitStreamArray(:); % Flatten the array

            % Ensure the bit stream array has the correct number of bits
            if length(bitStreamArray) > numBitsOriginal
                bitStreamArray = bitStreamArray(1:numBitsOriginal); % Trim extra bits
            elseif length(bitStreamArray) < numBitsOriginal
                error('Demodulated bit stream is shorter than expected. Something went wrong.');
            end

            % Assign to the output bit stream
            bit_stream(:, tx, disIdx) = bitStreamArray;

        end
    end
end
