% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

clear;
% Configuration Values
conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'

conf.f_s              = 48000;   % sampling rate  
conf.f_sym            = 100;     % symbol rate
conf.nframes          = 1;       % number of frames to transmit
conf.modulation_order = 2;       % BPSK:1, QPSK:2
conf.f_c              = 8000;

% OFDM Params
conf.spacing          = 9.375;
conf.nSubCarrier      = 256;
conf.ofdm_os_factor   = conf.f_s / (conf.spacing * conf.nSubCarrier);


conf.npreamble  = 100; 
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;

conf.training   = genpreamble(conf.nSubCarrier*conf.modulation_order);

% Init Section
% all calculations that you only have to do once
conf.os_factor  = conf.f_s/conf.f_sym;
if mod(conf.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end


% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.

guard_factors = 0.5;
n_OFDM_symbols = 300;
n_training_symbols = 1:1:5;

% Results
BER_cyclic = zeros(length(guard_factors), length(n_OFDM_symbols));
BER_training = zeros(length(n_training_symbols), length(n_OFDM_symbols));
BER_per_ofdm_diff_prefix = zeros(n_OFDM_symbols, length(guard_factors));
BER_per_ofdm_diff_num_train = zeros(n_OFDM_symbols, length(n_training_symbols));

% for k=1:conf.nframes
for n = 1:length(n_OFDM_symbols)
    fprintf('%d OFDM, ', n_OFDM_symbols(n));
    nOfdmSyms = n_OFDM_symbols(n);
    
    for j = 1:length(n_training_symbols)
        fprintf('%d training syms, ', n_training_symbols(j));
        nTrainSyms = n_training_symbols(j);
        conf.nTrainSyms = nTrainSyms; 
        % Adjust training symbols insertion indices
        nTotalOfdmSyms = nOfdmSyms + nTrainSyms;
        
        if nTrainSyms > 1
            conf.trainSymIdx = round(linspace(1,nTotalOfdmSyms,nTrainSyms));
        else
            conf.trainSymIdx = 1;
        end
        
        % Correct number of training symbols if provided input be arranged
        % conf.nTrainSyms = length(conf.trainSymIdx);
        
        conf.nOfdmSyms  = nOfdmSyms;
        conf.nbits      = conf.modulation_order * nOfdmSyms * conf.nSubCarrier;    % number of bits 
        conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

        for k = 1:length(guard_factors)
            fprintf('%f guard\n', guard_factors(k));
            % close all;
            % Use a different guard interval factor
            conf.guard_factor = guard_factors(k);

            % Generate random data
            txbits = randi([0 1],conf.nbits,1);

            % TODO: Implement tx() Transmit Function
            [txsignal, conf] = tx(txbits, conf, k);

            % % % % % % % % % % % %
            % Begin
            % Audio Transmission
            %

            % normalize values
            peakvalue       = max(abs(txsignal));
            normtxsignal    = txsignal / (peakvalue + 0.3);

            % create vector for transmission
            rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
            rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
            txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal

        %     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
            audiowrite('out.wav',rawtxsignal,conf.f_s)  

            % Platform native audio mode 
            if strcmp(conf.audiosystem,'native')

                % Windows WAV mode 
                if ispc()
                    disp('Windows WAV');
                    wavplay(rawtxsignal,conf.f_s,'async');
                    disp('Recording in Progress');
                    rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
                    disp('Recording complete')
                    rxsignal = rawrxsignal(1:end,1);

                % ALSA WAV mode 
                elseif isunix()
                    disp('Linux ALSA');
                    cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
                    system(cmd); 
                    disp('Recording in Progress');
                    system('aplay  out.wav')
                    pause(2);
                    disp('Recording complete')
                    rawrxsignal = wavread('in.wav');
                    rxsignal    = rawrxsignal(1:end,1);
                end

            % MATLAB audio mode
            elseif strcmp(conf.audiosystem,'matlab')
                disp('MATLAB generic');
                playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
                recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
                record(recobj);
                disp('Recording in Progress');
                playblocking(playobj)
                pause(0.5);
                stop(recobj);
                disp('Recording complete')
                rawrxsignal  = getaudiodata(recobj,'int16');
                rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;

            elseif strcmp(conf.audiosystem,'bypass')
                rawrxsignal = rawtxsignal(:,1);
                rxsignal    = rawrxsignal;
            end

            % Plot received signal for debgging
            figure;
            plot(rxsignal);
            title('Received Signal')

            %
            % End
            % Audio Transmission   
            % % % % % % % % % % % %

            % TODO: Implement rx() Receive Function
            [rxbits conf]       = rx(rxsignal,conf);

            res.rxnbits(1)      = length(rxbits);   % k instead of 1
            res.biterrors(1)    = sum(rxbits ~= txbits);   % k instead of 1
            ber = sum(res.biterrors)/sum(res.rxnbits)
            BER_cyclic(k,n) = ber;

            for i = 1:conf.nOfdmSyms
                r = rxbits((i-1)*conf.modulation_order*conf.nSubCarrier+1:(i)*conf.modulation_order*conf.nSubCarrier);
                t = txbits((i-1)*conf.modulation_order*conf.nSubCarrier+1:(i)*conf.modulation_order*conf.nSubCarrier);
                b = sum(r~=t)/length(r);
                fprintf('BER for OFDM symbol #%d: %.2f%%\n', i, b*100);
                BER_per_ofdm_diff_prefix(i,k)=b;
                BER_per_ofdm_diff_num_train(i,j)=b;
            end
        end
    end
end

%% %%%% multiple number of training symbols
x = 1:nOfdmSyms;
for i = 1:length(n_training_symbols)
    plot(x, BER_per_ofdm_diff_num_train(:,i) * 100);
    hold on;
end
hold off;
ylabel('BER (%)')
xlabel('nth OFDM Data Symbol')
title(sprintf('BER Plots for Each OFDM Symbol for Various Number of\nTraining Symbols with Fixed Cyclic Prefix Length 0.5 (Noisy Background)'))

for n = 1:length(n_training_symbols)
    leg{n} = sprintf('%d Training Symbols', n_training_symbols(n));
end
legend(leg);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%% 1 training symbol, multiple cyclic prefix lengths
x = 1:nOfdmSyms;
for i = 1:length(guard_factors)
    plot(x, BER_per_ofdm_diff_prefix(:,i) * 100);
    hold on;
end
hold off;
ylabel('BER (%)')
xlabel('nth OFDM Symbol After a Training Symbol')
title(sprintf('BER Plots for Each OFDM Symbol\nVarying Cyclic Prefix Length (Noisy Background)'))

for n = 1:length(guard_factors)
    leg{n} = sprintf('Cyclic Prefix Factor %.1f', guard_factors(n));
end
legend(leg);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%
plot(n_OFDM_symbols, BER_cyclic * 100);
ylabel('BER (%)')
xlabel('Number of OFDM Symbols')
title(sprintf('BER Plot for Various Number of OFDM Symbols\nwith Fixed Cyclic Prefix Length'))

for n = 1:length(guard_factors)
    leg{n} = sprintf('Cyclic Prefix Factor %.1f', guard_factors(n));
end
legend(leg);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%
plot(guard_factors, BER_cyclic * 100);
xlabel('Cyclic Prefix Length Factor')
title('BER Plot for Various Cyclic Prefix Lengths')

for n = 1:length(n_OFDM_symbols)
    leg{n} = sprintf('%d OFDM syms', n_OFDM_symbols(n));
end
legend(leg);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%
%per = sum(res.biterrors > 0)/conf.nframes
%ber = sum(res.biterrors)/sum(res.rxnbits)
