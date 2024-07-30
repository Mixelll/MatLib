function roundedValue = roundToPrecision(value, precision)
    % Rounds 'value' to 'precision' decimal places
    factor = 10^precision;
    roundedValue = round(value * factor) / factor;
end
