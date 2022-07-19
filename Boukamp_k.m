function[res] = Boukamp_k(r_dis_and_delay)
% Fit of IE-GPE data based on Boukamp's paper (1994).
% This function returns the sum of the squared residuals between the
% model values, evaluated using r_dis, the data and other parameters.  
% The function is called with a line such as:
% r_dis_fit = fminsearch('Boukamp_k',r_dis_guess);
% NB: no delay is fitted...
% many global variables are used:
global frac_18_O2_inf norm_36_t0 p corr_time q r_s S A B norm_36 norm_36_calc low_index high_index delta_f18; %frac_18_O2_t0
r_dis=r_dis_and_delay(1);
delay=r_dis_and_delay(2);
p=r_dis*S/A;
q=r_s*S*(1/A+1/B);
norm_36_calc=frac_18_O2_inf^2+(norm_36_t0-frac_18_O2_inf^2).*exp(-p*(corr_time-delay))-2*frac_18_O2_inf*(delta_f18).*(exp(-p*(corr_time-delay))-exp(-q*(corr_time-delay)))-(1-r_s/r_dis*(1+A/B))^2/(1-2*r_s/r_dis*(1+A/B))*(delta_f18)^2.*(exp(-p*(corr_time-delay))-exp(-q*(corr_time-delay)));
res = minus(norm_36(low_index:high_index,1),norm_36_calc(low_index:high_index,1));
res = sum(res.^2);
end