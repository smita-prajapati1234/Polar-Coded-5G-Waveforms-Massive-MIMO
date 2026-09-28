function encoded_bits = polar_encode(input_bits, N)
    % Polar Encoding
    % INPUTS:
    % input_bits: Input bit sequence as a row vector of 0s and 1s
    % N: Length of the encoded bit sequence (must be a power of 2)
    %
    % OUTPUT:
    % encoded_bits: Encoded bit sequence as a row vector of 0s and 1s

    % Check if N is a power of 2
    if bitand(N, N - 1) ~= 0
        error('N must be a power of 2 for polar encoding.');
    end

    % Calculate the number of bits to be frozen
    K = N / 2;
    
    % Calculate the number of information bits
    E = length(input_bits);

    % Extend the input sequence to the required length
    u = zeros(1, N);
    u(1:E) = input_bits;

    % Perform the polar encoding algorithm
    while K >= 1
        for i = 1:K
            u(i) = mod(u(i) + u(i + K), 2);
        end
        K = K / 2;
    end

    % The encoded bits are in the first E positions of the vector u
    encoded_bits = u(1:E);
end
