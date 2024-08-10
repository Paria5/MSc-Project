function modulated_signal = modulation(bit_stream, modulation_type, M, codingRate)
    [numBits, numTx] = size(bit_stream);
    bitsPerSymbol = log2(M);
    
    % Validate inputs
    validateattributes(modulation_type, {'char', 'string'}, {'nonempty'});
    validateattributes(M, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(codingRate, {'numeric'}, {'scalar', '>', 0, '<=', 1});

    % Apply coding to bit stream
    coded_bits = [];
    for tx = 1:numTx
        encoded_stream = encode_bitstream(bit_stream(:, tx).', codingRate); % Encode each transmitter's bit stream
        coded_bits = [coded_bits; encoded_stream];
    end

   % Ensure the coded bit stream is a row vector
    coded_bits = coded_bits(:).'; % Flatten and ensure it's a row vector

    % Initialize matrix for modulated signals
    numCodedBits = length(coded_bits);
    numSymbols = ceil(numCodedBits / bitsPerSymbol);
    padded_bs = [coded_bits, zeros(1, numSymbols * bitsPerSymbol - numCodedBits)]; % Row vector padding
    modulated_signal = zeros(numSymbols, numTx); 

    % Process each column of the bit stream
    for tx = 1:numTx
        % Extract the current column (transmitter)
        current_bit_stream = padded_bs;
        reshape_bs = reshape(current_bit_stream, bitsPerSymbol, []);
        dataStream = bi2de(reshape_bs.', 'left-msb');

        % Perform modulation based on type
        if strcmp(modulation_type, 'psk') % PSK modulation
            modulated_signal(:, tx) = pskmod(dataStream, M);
        elseif strcmp(modulation_type, 'qam') % QAM modulation
            modulated_signal(:, tx) = qammod(dataStream, M);
        else
            error('Unsupported modulation type.');
        end

    end
end
