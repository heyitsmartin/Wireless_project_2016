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

preamble = lfsr_framesync(conf.npreamble);
signal = qpsk_mapper([preamble; txbits]); 
signal = upsample(signal, os_factor);

time = linspace(1, 1+length(signal)/f_s, length(signal));

sig_im = sin(2*pi*f_c*time);
sig_re = cos(2*pi*f_c*time);

ps_filter = rrc(os_factor, 0.22, 20);

txsignal = conv(signal, ps_filter, 'same');
txsignal = txsignal .* transp(exp(2*pi*1i*f_c*time));
% txsignal = sig_re'.*real(txsignal) - sig_im'.*imag(txsignal);




