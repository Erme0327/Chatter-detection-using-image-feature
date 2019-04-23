% data processing -to do FFT, using frequency averaging, 
% This code was developed by Huaizhong Li in Griffith University, Australia

 
function [data_fft, f_delta_act] = fft_yun(data_input, f_delta, f_max, f_sampling)
% the following parameters are assumed:
% data_input: input data to do fft;
% f_delta: expected frequency resolution but not the final one. 
% the actual freq resolution will be returned by f_delta_act, 
% f_max: expected maximum frequency (range), 
% f_sampling: sampling frequency of input data in Hz.
% data_fft: return the fft results.

Fs = f_sampling;
n_total = length(data_input); 
 
%calculate the data size and actural frequency:
T_needed = ceil(Fs*1/f_delta);
n_length = pow2(nextpow2(T_needed));  % Transform length
% optimising speed by considering segment size as 2^n  

if n_length > n_total
    n_length = n_total;
end
 
f_delta_act = Fs / n_length;
n_f_max = round(f_max / f_delta_act);
if n_f_max >n_length/2
    n_f_max = n_length/2;
end
 
%create a Hanning window of n_length points:
w=hann(n_length);

% devide the raw data into a few segments, then using energy averaging
% note: overlapping 0
% n_segment = floor(n_total/n_length);
ratio_overlap = 0; % setting the overlap ratio, 0-1.
n_segment = 1+ floor((n_total-n_length )/((1-ratio_overlap)*n_length));
%  n_length,n_total,n_segment
pos = 1;
% n_order = 1;
 
for n_order = 1:n_segment 
    tmp = data_input(pos:pos+n_length-1);
    x_data =tmp-mean(tmp);
    % apply the window:
    x_data = x_data.*w;
    y = fft(x_data,n_length);           % DFT
    data_amp(:,n_order) = abs(y(1:n_f_max)/(n_length))*2;% two sided FFT, so multiply by 2 
    
    pos = pos+floor((1-ratio_overlap)*n_length);
end
data_fft = mean(data_amp,2);
 
