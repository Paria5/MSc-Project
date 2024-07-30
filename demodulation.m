function bit_stream = demodulation(modulated_signal, modulation_type, M)
    [numSymbols, numTx, numUsers] = size(modulated_signal);
    bitsPerSymbol = log2(M);
    
    % Validate inputs
    validateattributes(modulation_type, {'char', 'string'}, {'nonempty'});
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'positive'});
    
    % Initialize matrix for demodulated bits
    numBits = numSymbols * bitsPerSymbol;
    bit_stream = zeros(numBits, numTx, numUsers); 

    % Process each user's signal
    for userIdx = 1:numUsers
        for tx = 1:numTx
            % Extract the current column (transmitter) for the current user
            current_modulated_signal = modulated_signal(:, tx, userIdx);

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
            if length(bitStreamArray) > numBits
                bitStreamArray = bitStreamArray(1:numBits); % Trim extra bits
            elseif length(bitStreamArray) < numBits
                bitStreamArray = [bitStreamArray; zeros(numBits - length(bitStreamArray), 1)]; % Pad with zeros if necessary
            end

            % Assign to the output bit stream
            bit_stream(:, tx, userIdx) = bitStreamArray;

            % Debugging information
            disp(['User ', num2str(userIdx), ', Tx ', num2str(tx), ':']);
            disp(['Expected number of bits: ', num2str(numBits)]);
            disp(['Actual number of bits: ', num2str(length(bitStreamArray))]);
        end
    end
end
