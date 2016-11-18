function symbols = map(bits)

% Map symbols
symbols = ((2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1));

% Normalize power
symbols = symbols/(sqrt(2));