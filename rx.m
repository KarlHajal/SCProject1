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

    f_c = conf.f_c;
    f_s = conf.f_s;
    nbits = conf.nbits;
    modulation_order = conf.modulation_order;
    number_of_symbols = nbits / modulation_order;
    os_factor = conf.os_factor;
    npreamble = conf.npreamble;

    % dummy 
    % rxbits = zeros(conf.nbits,1);

    % down conversion / de-multiplexing
    carrier_seq = exp(-1i*2*pi*(f_c/f_s)*(1:length(rxsignal))).';
    rxsyms_dc = rxsignal.*carrier_seq;
    % low pass filter
    rxsyms_lp = 2*lowpass(rxsyms_dc, conf);
    % matched filter
    filtered_rxsyms = conv(rxsyms_lp, rrc(os_factor), 'same');
    % frame synchronizationg
    [data_idx, theta, magnitude] = frame_sync(filtered_rxsyms, conf);

    
    cum_err = 0;
    % Use preamble symbols to improve timing offset estimation
    for i = 1 : npreamble
        idx_start = data_idx - npreamble*os_factor + (i-1)*os_factor;
        idx_range = idx_start : idx_start+os_factor-1;
        segment = filtered_rxsyms(idx_range);
        [~, cum_err] = epsilon_estimation(segment, cum_err);
    end
    
    theta_hat = theta;
    data = zeros(1, number_of_symbols);
    for i = 1 : number_of_symbols
         idx_start  = data_idx + (i-1)*os_factor;
         idx_range  = idx_start : idx_start+os_factor-1;
         segment    = filtered_rxsyms(idx_range);
         
         [epsilon, cum_err] = epsilon_estimation(segment, cum_err);
         
         data(i) = cubic_interpolation(epsilon, filtered_rxsyms, os_factor, idx_start);
         %data(i) = linear_interpolation(epsilon, filtered_rxsyms, os_factor, idx_start);
         
         % Phase estimation    
         % Apply viterbi-viterbi algorithm
         deltaTheta = 1/4*angle(-data(i)^4) + pi/2*(-1:4);
         
         [~, ind] = min(abs(deltaTheta - theta_hat));
         theta = deltaTheta(ind);
         theta_hat = mod(0.01*theta + 0.99*theta_hat, 2*pi);
         data(i) = (1/magnitude) * data(i) * exp(-1j * theta_hat);
    end
    
    rxbits = demapper(data);
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
    %pwr = abs(segment).^2;
    %tmp_seq = (-1i).^(0 : length(segment)-1);
    %diff_err = tmp_seq*pwr;
    result = fft(abs(segment).^2);
    diff_err = result(2);
    cum_err = cum_err + diff_err;
    epsilon = -1/(2*pi)*angle(cum_err);
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