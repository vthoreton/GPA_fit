function[res] = Den_Otter_two_steps(input_guess)
% This function returns the sum of the squared residuals between the
% normalised data (32, 36, 34 weighted) and the model values,
%evaluated in the range (low:high), using p.  
% The function is called with a line such as:
% p_calc=fminsearch('fit_two_steps_model',input_guess);
% two other global variables, two_steps_model_data and fitting_param must be define and contain: {corr_time,norm_32,norm_34,norm_36,A,B,low_DO,high_DO,tau_calc};
%p can be either is either a single number (p1=p2) or a double (p1 and p2)
% it can also be a triplet then tau (3rd parameter) is refined too
global two_steps_model_data fitting_param frac_18_O2_t0 frac_18_O2_inf f18_sample_t0 norm_36_t0 norm_36_inf
%ponderation coefficient for mass 34 is set to 1:
pond_34=1;
time=two_steps_model_data{1};
norm_32=two_steps_model_data{2};
norm_34=two_steps_model_data{3};
norm_36=two_steps_model_data{4};
A=two_steps_model_data{5};
B=two_steps_model_data{6};
low=two_steps_model_data{7};
high=two_steps_model_data{8};

p1_set=fitting_param(1);
p2_set=fitting_param(2);
delay_set=fitting_param(3);
tau_set=fitting_param(4);
f_mod=fitting_param(5);

n=1;
if bitget(f_mod,4)==0
p1=p1_set;
elseif bitget(f_mod,4)==1 & bitget(f_mod,3)==0
p1=input_guess(n);
p2=p1;
n=n+1;
else
p1=input_guess(n);
n=n+1;  
end
if bitget(f_mod,3)==0 & bitget(f_mod,4)==0
p2=p2_set;
elseif bitget(f_mod,3)==0 & bitget(f_mod,4)==1
%case already considered
else
p2=input_guess(n);
n=n+1;
end
if bitget(f_mod,2)==0
delay=delay_set;
else
delay=input_guess(n);
n=n+1;
end
if bitget(f_mod,1)==0
tau_1=tau_set;
else
tau_1=input_guess(n);
end

        f_36_t0=norm_36_t0;
        f_36_inf=norm_36_inf;
        f_18g_t0=frac_18_O2_t0;
        f_18b_t0=f18_sample_t0;
        f_18_inf=frac_18_O2_inf;
        beta_0=f_18_inf^2;
        beta_1=f_18_inf*((f_18b_t0-f_18g_t0)*(p1+p2)+2*(f_18g_t0-f_18_inf));
        beta_2=((f_18b_t0-f_18g_t0)*p1+f_18g_t0-f_18_inf)*((f_18b_t0-f_18g_t0)*p2+f_18g_t0-f_18_inf);
        tau_2=tau_1*(1-beta_1/(2*f_18_inf*(f_18g_t0-f_18_inf)));
        epsilon_1=2*f_18_inf*(f_18g_t0-f_18_inf);
        epsilon_2=beta_2*tau_1/(tau_1-2*tau_2);
        f_36=f_18_inf^2+epsilon_1*exp(-(time-delay)/tau_1)+epsilon_2*exp(-2*(time-delay)/tau_1)+(f_36_t0-f_18_inf^2-epsilon_1-epsilon_2)*exp(-(time-delay)/tau_2);
        f_18g=f_18_inf+(f_18g_t0-f_18_inf)*exp(-(time-delay)/tau_1);
        f_18b=f_18_inf+(f_18b_t0-f_18_inf)*exp(-(time-delay)/tau_1);
        f_34=2*(f_18g-f_36);
        f_32=1-f_34-f_36;
        %res =minus(f_34(low:high,1),norm_34(low:high,1));
        res = sum((minus(f_36(low:high,1),norm_36(low:high,1))).^2)+pond_34*sum((minus(f_34(low:high,1),norm_34(low:high,1))).^2)+sum((minus(f_32(low:high,1),norm_32(low:high,1))).^2);
end