function bit_stream = demodulation(modulated_signal, modulation_type, M,numBitsOriginal, codingRate)
    [numSymbols, numTx, numDistance] = size(modulated_signal);
    bitsPerSymbol = log2(M);
    
    % Validate inputs
    validateattributes(modulation_type, {'char', 'string'}, {'nonempty'});
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(numBitsOriginal, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(codingRate, {'numeric'}, {'scalar', '>', 0, '<=', 1});

    % Calculate the repeat factor based on the coding rate
    repeatFactor = round(1 / codingRate);

    % Initialize matrix for demodulated bits
    bit_stream = zeros(numBitsOriginal, numTx, numDistance); 

    % Process each user's signalh
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


            % Decoding: Reverse the repetition by majority voting or selecting the first bit
            decodedBits = zeros(numBitsOriginal, 1); % Initialize decoded bit stream
            for i = 1:numBitsOriginal
                repeatedBits = bitStreamArray((i-1)*repeatFactor+1:i*repeatFactor);
                decodedBits(i) = mode(repeatedBits); % Majority voting (can also take first bit)
            end

            % Ensure the decoded bit stream has the correct number of bits
            if length(decodedBits) > numBitsOriginal
                decodedBits = decodedBits(1:numBitsOriginal); % Trim extra bits
            elseif length(decodedBits) < numBitsOriginal
                error('Decoded bit stream is shorter than expected. Something went wrong.');
            end

            % Assign to the output bit stream
            bit_stream(:, tx, disIdx) = decodedBits;

        end
    end
end
