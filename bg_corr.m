function [res] = bg_corr(bg_guess)
global bg_corr_data 

bg_intensity32=bg_guess(1);
bg_intensity34=bg_guess(2);
bg_intensity36=bg_guess(3);
time=bg_corr_data{1};
raw_intensity32=bg_corr_data{2};
raw_intensity34=bg_corr_data{3};
raw_intensity36=bg_corr_data{4};
frac_18_O2_t0=bg_corr_data{5};
frac_18_O2_inf=bg_corr_data{6};
tau=bg_corr_data{7};
delay=bg_corr_data{8};
low=bg_corr_data{9};
high=bg_corr_data{10};



corr_intensity32=raw_intensity32-bg_intensity32;
corr_intensity34=raw_intensity34-bg_intensity34;
corr_intensity36=raw_intensity36-bg_intensity36;

sum_intensityO2=corr_intensity32+corr_intensity34+corr_intensity36;
norm_32=corr_intensity32./sum_intensityO2;
norm_34=corr_intensity34./sum_intensityO2;
norm_36=corr_intensity36./sum_intensityO2;
frac_18_O2=norm_36+1/2*norm_34;
norm_frac_18_O2=(frac_18_O2-frac_18_O2_t0)/(frac_18_O2_inf-frac_18_O2_t0);

norm_calc=1-exp(-1/tau.*(time-delay));

res = minus(norm_frac_18_O2(low:high,1),norm_calc(low:high,1));
res = sum(res.^2);
end