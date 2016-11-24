function [rxbits conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%
%phase correction
f_c = conf.f_c;
f_s = conf.f_s;
nsyms = conf.nsyms;
os_factor = conf.os_factor;
modulation_order = conf.modulation_order; % BPSK:1, QPSK:2

time = linspace(1, 1+length(rxsignal)/f_s, length(rxsignal));
rx_comp = rxsignal .* transp(exp(-2*1i*pi*f_c*time));       % down freq conv

rx_bb = 2 * lowpass(rx_comp, conf);

ps_filter = rrc(os_factor, 0.22, 20);
rx_matched = conv(rx_bb, ps_filter, 'same');
[beginning_of_data, phase_of_peak, ~] = frame_sync(rx_matched, os_factor);

[rxbits, ~] = phase_correction(rx_matched, nsyms, beginning_of_data, os_factor, phase_of_peak);

if (modulation_order == 1)
    rxbits = bpsk_demapper(rxbits);    
elseif (modulation_order == 2)
    rxbits = qpsk_demapper(rxbits);    
end


