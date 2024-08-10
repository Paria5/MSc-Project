function encoded_stream = encode_bitstream(bit_stream, codingRate)
    % Example simple coding: repeat each bit 1/codingRate times
    repeatFactor = round(1 / codingRate);
    encoded_stream = repelem(bit_stream, repeatFactor);
end
