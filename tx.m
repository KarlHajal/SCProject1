function [txsignal conf] = tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

% dummy 400Hz sinus generation
% time = 1:1/conf.f_s:4;
% txsignal = 0.3*sin(2*pi*400 * time.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modulation_order = conf.modulation_order;
os_factor = conf.os_factor;
f_c = conf.f_c;
f_s = conf.f_s;
npreamble = conf.npreamble;

preamble_sym = modulator(preamble_generate(npreamble), 1);

txsyms = [preamble_sym; modulator(txbits, modulation_order)]; % bits to symbols;
% scatter(real(txsyms), imag(txsyms), 'filled');
txsyms_us = upsample(txsyms, os_factor); % upsampling
txsyms_rrc = conv(txsyms_us, rrc(os_factor), 'same'); % pulse shaping
carrier_seq = exp(1i*2*pi*(f_c/f_s)*(1:length(txsyms_rrc))).'; % carrier sequence
txsyms_uc = txsyms_rrc.*carrier_seq; % up conversion
txsignal = real(txsyms_uc);
figure(2); plot(txsignal(npreamble*os_factor+1:(npreamble+5.5)*os_factor), 'Linewidth', 1); ylabel('txsignal')

end

