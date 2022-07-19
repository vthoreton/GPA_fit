function [res] = exponential_decay_extract_tau_and_delay(input_guess)
% This function returns the sum of the squared residuals between the
% normalised data (NB:adjusted on 0) (time, norm_data) and the model values,
%evaluated in the range (low:high), using tau_and_delay_guess that contains
%estimate values of both tau and a delay.
% The function is called with a line such as:
% tau_calc = fminsearch('extract_tau_and_delay',guess);
% An global variable, extract_tau_data must be defined and contains time,norm_data,low,high
global extract_tau_data fitting_param;
time=extract_tau_data{1};
norm_data=extract_tau_data{2};
low=extract_tau_data{3};
high=extract_tau_data{4};

tau_set=fitting_param(1);
delay_set=fitting_param(2);
f_mod=fitting_param(3);
n=1;
if bitget(f_mod,2)==0
tau=tau_set;
else
tau=input_guess(n);
n=n+1;
end
if bitget(f_mod,1)==0
delay=delay_set;
else
delay=input_guess(n);
end

norm_calc=1-exp(-1/tau.*(time-delay));
res = minus(norm_data(low:high,1),norm_calc(low:high,1));
res = sum(res.^2);
end