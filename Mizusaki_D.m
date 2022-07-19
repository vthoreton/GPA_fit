function [res] = Mizusaki_D(input_guess)
% Mizusaki_D: fitting D (as well as a delay and possibly the particle radius) for GPA
% requires access to one global variable:
%   - fitting_data containing the data [time norm18 A B radius low high]
% fitting mode is defined by the number of parameters in the input:
% 1*fitting D only
% 2*fitting D and a delay
% 3*fitting D,a delay and the particle size !Do not use this last mode
% especially if you plan to use Mizusaki_kD after!
syms n x;
global fitting_data;

t=fitting_data{1};
norm18=fitting_data{2};
A=fitting_data{3};
B=fitting_data{4};
radius=fitting_data{5};%sphere radius (m)
low=fitting_data{6};
high=fitting_data{7};
q=fitting_data{8};
alpha=A/B;
M=zeros(size(t));

if size(input_guess,2)==1
    D=input_guess;
    delay=0;
elseif size(input_guess,2)==2 
    D=input_guess(1);
    delay=input_guess(2);
elseif size(input_guess,2)==3 
    D=input_guess(1);
    delay=input_guess(2);
    radius=input_guess(3);    
end
%sum:
for n=1:1:size(q,1)
M=M+(6*alpha*(alpha+1)*exp(-D*q(n)^2.*(t-delay)/radius^2))/(9+9*alpha+q(n)^2*alpha^2);
end
M=1-M;
res = minus(M(low:high,1),norm18(low:high,1));
res = sum(res.^2);
end