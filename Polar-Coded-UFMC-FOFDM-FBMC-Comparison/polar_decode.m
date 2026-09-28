
function decoded_bits = polar_decode(received_bits, N)
    % Polar Decoding
    % INPUTS:
    % received_bits: Received bit sequence as a row vector of 0s and 1s
    % N: Length of the received bit sequence (must be a power of 2)
    %
    % OUTPUT:
    % decoded_bits: Decoded bit sequence as a row vector of 0s and 1s

    % Check if N is a power of 2
    if bitand(N, N - 1) ~= 0
        error('N must be a power of 2 for polar decoding.');
    end

    % Calculate the number of bits to be frozen
    K = N / 2;
    
    % Calculate the number of information bits
    E = length(received_bits);

    % Extend the received sequence to the required length
    y = zeros(1, N);
    y(1:E) = received_bits;

    % Perform the polar decoding algorithm
    while K >= 1
        for i = 1:K
            y(i) = mod(y(i) + y(i + K), 2);
        end
        K = K / 2;
    end

    % The decoded bits are in the first E positions of the vector y
    decoded_bits = y(1:E);
end
