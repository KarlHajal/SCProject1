function out_syms = modulator(in_bits, modulation_order)

if (modulation_order~=1 && modulation_order~=2)
    error('Modulation_order wrong!')
end

if modulation_order == 1 % BPSK
    out_syms = -2*in_bits+1; % 0 -> -1; 1 -> 1 
else % QPSK
    out_syms = zeros(ceil(size(in_bits, 1)/modulation_order), 1);
    GrayMap = 1/sqrt(2)*[-1-1i, -1+1i, 1-1i, 1+1i];
    for i = 1:size(out_syms, 1)
        decimal_index = 2*in_bits(2*i-1)+in_bits(2*i);
        out_syms(i) = GrayMap( decimal_index + 1 );
    end
end

end