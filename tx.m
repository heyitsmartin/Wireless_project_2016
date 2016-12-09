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

f_c = conf.f_c;
f_s = conf.f_s;
os_factor = conf.os_factor;
ofdm_os_factor = conf.ofdm_os_factor;
modulation_order = conf.modulation_order; % BPSK:1, QPSK:2
nOfdmSyms = conf.nOfdmSyms;
nSubCarrier = conf.nSubCarrier;
ofdm_signal = [];

% User our lfsr to generate a training OFDM symbol.
txbits = [conf.training; txbits];

preamble = 1 - 2*genpreamble(conf.npreamble); % BPSK

if (modulation_order == 1)
    mapped = bpsk_mapper(txbits);
elseif (modulation_order == 2)
    mapped = qpsk_mapper(txbits);
end

ofdm_start_idx = 1;
for i = 1:nOfdmSyms+1
    mapped_symbols = mapped(ofdm_start_idx:ofdm_start_idx+nSubCarrier-1);
    ofdm_symbol = osifft(mapped_symbols, ofdm_os_factor);
    ofdm_signal = [ofdm_signal; ofdm_symbol(1/2*length(ofdm_symbol)+1:end); ofdm_symbol]; % ofdm symbol with guard interval
    ofdm_start_idx = ofdm_start_idx + nSubCarrier;
end

% signal = [preamble; mapped];       % bpsk preamble
% signal = upsample(signal, os_factor);

ps_filter = rrc(os_factor, 0.22, 10*os_factor);
% upsample preamble separately because we don't perform IDFT on it.
preamble_ = upsample(preamble, os_factor); 
preamble = conv(preamble_, ps_filter, 'same');

signal = [preamble; ofdm_signal];
% txsignal = conv(signal, ps_filter, 'same');

time = linspace(1, 1+length(signal)/f_s, length(signal));

sig_re = cos(2*pi*f_c*time);
sig_im = sin(2*pi*f_c*time);

txsignal = sig_re'.*real(signal) - sig_im'.*imag(signal);
