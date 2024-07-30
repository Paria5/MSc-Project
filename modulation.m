function modulated_signal = modulation(bit_stream, modulation_type, M)
    [numBits, numTx] = size(bit_stream);
    bitsPerSymbol = log2(M);
    
    % Validate inputs
    validateattributes(modulation_type, {'char', 'string'}, {'nonempty'});
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'positive'});

    % Initialize matrix for modulated signals
    numSymbols = ceil(numBits / bitsPerSymbol);
    modulated_signal = zeros(numSymbols, numTx); 

    % Process each column of the bit stream
    for tx = 1:numTx
        % Extract the current column (transmitter)
        current_bit_stream = bit_stream(:, tx).'; % Transpose to make it a row vector

        % Perform modulation based on type
        if strcmp(modulation_type, 'psk') % PSK modulation
            % Pad bit stream to ensure it fits an integer number of symbols
            padded_bs = [current_bit_stream, zeros(1, numSymbols * bitsPerSymbol - numBits)]; % Row vector padding
            reshape_bs = reshape(padded_bs, bitsPerSymbol, []);
            dataStream = bi2de(reshape_bs.', 'left-msb');
            modulated_signal(:, tx) = pskmod(dataStream, M);

        elseif strcmp(modulation_type, 'qam') % QAM modulation
            % Pad bit stream to ensure it fits an integer number of symbols
            padded_bs = [current_bit_stream, zeros(1, numSymbols * bitsPerSymbol - numBits)]; % Row vector padding
            reshape_bs = reshape(padded_bs, bitsPerSymbol, []);
            dataStream = bi2de(reshape_bs.', 'left-msb');
            modulated_signal(:, tx) = qammod(dataStream, M);

        else
            error('Unsupported modulation type.');
        end
    end
end