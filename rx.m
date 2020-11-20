function [rxbits conf] = rx(rxsignal, txsyms, conf,k)
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

    f_c = conf.f_c;
    f_s = conf.f_s;
    nbits = conf.nbits;
    modulation_order = conf.modulation_order;
    number_of_symbols = ceil(nbits / modulation_order);
    os_factor = conf.os_factor;
    npreamble = conf.npreamble;

    % dummy 
    % rxbits = zeros(conf.nbits,1);

    filtered_rxsyms = preprocess(rxsignal, conf);
    
    % frame synchronization
    [data_idx, theta, magnitude] = frame_sync(filtered_rxsyms, conf);

    
    cum_err = 0;
    % Use preamble symbols to improve timing offset estimation
    for i = 1 : npreamble
        idx_start = data_idx - npreamble*os_factor + (i-1)*os_factor;
        idx_range = idx_start : idx_start+os_factor-1;
        segment = filtered_rxsyms(idx_range);
        [~, cum_err] = epsilon_estimation(segment, cum_err);
    end
    
    data_time = zeros(ceil(nbits/modulation_order), 1);
    data_phase = zeros(ceil(nbits/modulation_order), 1);
    epsilon_arr = zeros(ceil(nbits/modulation_order), 1);
    theta_hat = zeros(ceil(nbits/modulation_order)+1, 1); 
    theta_hat(1) = theta;    
    
    %theta_hat = theta;
    %data = zeros(1, number_of_symbols);
    for i = 1 : number_of_symbols
         idx_start  = data_idx + (i-1)*os_factor;
         idx_range  = idx_start : idx_start+os_factor-1;
         segment    = filtered_rxsyms(idx_range);
         
         [epsilon_arr(i), cum_err] = epsilon_estimation(segment, cum_err);
         
         data_time(i) = cubic_interpolation(epsilon_arr(i), filtered_rxsyms, os_factor, idx_start);
         %data_time(i) = linear_interpolation(epsilon, filtered_rxsyms, os_factor, idx_start);
         
         % Phase estimation    
         % Apply viterbi-viterbi algorithm
         deltaTheta = 1/4*angle(-data_time(i)^4) + pi/2*(-1:4);
         
         [~, ind] = min(abs(deltaTheta - theta_hat(i)));
         theta = deltaTheta(ind);
         
         theta_hat(i+1) = mod(0.01*theta + 0.99*theta_hat(i), 2*pi);
         
         data_phase(i) =  data_time(i) * exp(-1j * theta_hat(i+1));
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the magnitude and phase of the estimated and original symbols
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %txsyms_mf = preprocess(txsyms, conf);
    %txsyms_sample = txsyms_mf(os_factor*npreamble+1:os_factor:length(txsyms_mf));
    %figure(3); 
    %subplot(2, 1, 1);
    %plot(abs(txsyms_sample(1:100)), 'Linewidth', 1.5);hold on; 
    %plot(abs(data_phase(1:100)), 'Linewidth', 1.5); hold off;
    %ylabel('magnitude'); legend('txsyms', 'rxsyms after correction');
    %subplot(2, 1, 2);
    %plot(angle(txsyms_sample(1:100)), 'Linewidth', 1.5); hold on; 
    %plot(angle(data_phase(1:100)), 'Linewidth', 1.5); hold off;
    %ylabel('phase'); legend('txsyms', 'rxsyms after correction');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    rxbits = demapper(data_phase);
end

function [out_syms] = preprocess(in_syms, conf)
    f_c = conf.f_c; f_s = conf.f_s; os_factor = conf.os_factor;
    carrier_seq = exp(-1i*2*pi*(f_c/f_s)*(1:length(in_syms))).';
    in_syms_dc = in_syms.*carrier_seq;
    in_syms_lp = 2*lowpass(in_syms_dc, conf);
    out_syms = conv(in_syms_lp, rrc(os_factor), 'same');
end


function [beginning_of_data, phase_of_peak, magnitude_of_peak] = frame_sync(in_syms, conf)

    npreamble = conf.npreamble;
    os_factor = conf.os_factor;

    preamble_syms = modulator(preamble_generate(npreamble), 1);
    current_peak_value = 0;
    samples_after_threshold = os_factor;
    detection_threshold = 15;

    for i = os_factor*npreamble+1:length(in_syms)
        r = in_syms(i-os_factor*npreamble : os_factor : i-os_factor); 
        c = preamble_syms'*r;
        T = abs(c)^2/abs((r')*r);
        
        if (T > detection_threshold || samples_after_threshold < os_factor)
            samples_after_threshold = samples_after_threshold - 1;
            if (T > current_peak_value)
                beginning_of_data = i;
                phase_of_peak = mod(angle(c),2*pi);
                magnitude_of_peak = abs(c)/npreamble;
                current_peak_value = T;
            end
            if (samples_after_threshold == 0)
                return;
            end
        end
    end
end

function [epsilon, cum_err] = epsilon_estimation(segment, cum_err)
    os_factor = length(segment);  
    tmp_seq = exp(-1i*(2*pi)*(0:os_factor-1)*(1/os_factor));
    cum_err = cum_err + tmp_seq*(abs(segment).^2); 
    epsilon = -(1/(2*pi))*angle(cum_err);
end

function [y_hat] = linear_interpolation(epsilon, filtered_rxsyms, os_factor, idx_start)
    sample_diff = floor(epsilon*os_factor); % integer
    int_diff = mod(epsilon*os_factor,1); % interval [0 1)
    y = filtered_rxsyms(idx_start+sample_diff : idx_start+sample_diff+1);
    y_hat = y(1)+int_diff*(y(2)-y(1));
end

function [y_hat] = cubic_interpolation(epsilon, filtered_rxsyms, os_factor, idx_start)
    sample_diff   = floor(epsilon*os_factor); % integer
    int_diff      = mod(epsilon*os_factor,1); % interval [0 1)
    y = filtered_rxsyms( idx_start+sample_diff-1 : idx_start+sample_diff+2 );
    c = 1/6 * [-1 3 -3 1 ; 3 -6 3 0 ; -2 -3 6 -1 ; 0 6 0 0] * y;
    y_hat = c(1) * int_diff^3 + c(2) * int_diff^2 + c(3) * int_diff + c(4);
end