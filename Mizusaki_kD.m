function [res] = Mizusaki_kD(input_guess)
% Mizusaki_kD: fitting of k and D (as well as a delay and the particle radius) for GPA
% requires access to two global variables:
%   - fitting_data containing the data [time norm18 A B radius low high]
%   ;radius is not used
%   - fitting_param containing the initial fitting parameters [D k delay radius f_mod]
% fitting_param contains the initial values of [D k delay radius] as well as
% a binary number of 4 bits f_mod telling which of D k delay and radius should be
% refined (1: refine, 0: keep constant)

% input_guess must contain the parameters among D k delay and radius that
% WILL be refined. It is important not to put any other (constant) parameter in.
%e.g. input_guess = [D k delay] comes with f_mod=0b1110
%e.g. input_guess = [k radius] comes with f_mod=0b0101 

%%%syms n x;
global fitting_data fitting_param q nroots;

% getting the data from fitting_data
t=fitting_data{1};
norm18=fitting_data{2};
A=fitting_data{3};
B=fitting_data{4};
low=fitting_data{6};
high=fitting_data{7};

% getting fitting parameters fitting_param
D_set=fitting_param(1);
k_set=fitting_param(2);
delay_set=fitting_param(3);
radius_set=fitting_param(4);
f_mod=fitting_param(5);


n=1;
if bitget(f_mod,4)==0
D=D_set;
else
D=input_guess(n);
n=n+1;
end
if bitget(f_mod,3)==0
k=k_set;
else
k=input_guess(n);
n=n+1;
end
if bitget(f_mod,2)==0
delay=delay_set;
else
delay=input_guess(n);
n=n+1;
end
if bitget(f_mod,1)==0
radius=radius_set;
else
radius=input_guess(n);
end

alpha=A/B;
M=zeros(size(t));
kaD=k*radius/D;
q=Mizusaki_kD_roots([k D alpha radius nroots]);
max=nroots;
if size(q,1) < nroots
   disp('Warning! Iteration using less roots than specified.')
   size(q,1)
   max=size(q,1);
end
for n=1:1:max
M=M+(6*alpha*kaD^2*(alpha+1)*exp(-D*q(n)^2.*(t-delay)/radius^2))/(alpha^2*q(n)^4+alpha*kaD*(alpha*(kaD-1)-6)*q(n)^2+9*(1+alpha)*kaD^2);
end
M=1-M;
res = minus(M(low:high,1),norm18(low:high,1));
res = sum(res.^2);
end