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
ofdm_os_factor = conf.ofdm_os_factor;
nOfdmSyms = conf.nOfdmSyms;
nSubCarrier = conf.nSubCarrier;
nTrainSyms = conf.nTrainSyms;
spacing = conf.spacing;
channel_offsets = zeros(nSubCarrier, 1);
rxbits = [];

time = linspace(1, 1+length(rxsignal)/f_s, length(rxsignal));
rx_comp = rxsignal .* transp(exp(-2*1i*pi*f_c*time));       % down freq conv

rx_bb = ofdmlowpass(rx_comp, conf, 2000);

ps_filter = rrc(os_factor, 0.22, 10*os_factor);
rx_matched = conv(rx_bb, ps_filter, 'same');
[beginning_of_data, ~, ~] = frame_sync(rx_matched, os_factor);
window_len = nSubCarrier * ofdm_os_factor;
for i = 1:nOfdmSyms+nTrainSyms
    beginning_of_data = beginning_of_data + floor(conf.guard_factor * window_len); % skip over cyclic prefix
    end_of_symbol = beginning_of_data + window_len - 1;
    ofdm_symbol = rx_bb(beginning_of_data:end_of_symbol);
    mapped_bits = osfft(ofdm_symbol, ofdm_os_factor);
    
    a = sum(conf.training ~= qpsk_demapper(mapped_bits));
    if sum(i == conf.trainSymIdx) > 0 
        % If trainSymIdx includes this OFDM symbol index, its a training symbol
        training = conf.training;
        channel_offsets = mapped_bits ./ qpsk_mapper(training);
    else
        mapped_bits = mapped_bits ./ channel_offsets;
        rxbits = [rxbits; mapped_bits];
        % Constellation plot
%         if (mod(i,50) == 2)
%             scatter(real(mapped_bits), imag(mapped_bits), 'MarkerFaceColor', rand(1,3));
%             figure
%             title(sprintf('i = %d',i));
%         end
    end
    beginning_of_data = end_of_symbol + 1;
end

% [rxbits, ~] = phase_correction(rx_matched, nsyms, beginning_of_data, os_factor, phase_of_peak);

if (modulation_order == 1)
    rxbits = bpsk_demapper(rxbits);
elseif (modulation_order == 2)
    rxbits = qpsk_demapper(rxbits);
end
