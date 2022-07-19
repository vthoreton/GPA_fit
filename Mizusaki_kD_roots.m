function[q] = Mizusaki_kD_roots(k_D_alpha_radius_nroots)
%Determine the roots of tan(x)=(3M*x-alpha*x^3)/(3M+(M-1)*alpha*x^2)
sym x;
global xroots;
k=k_D_alpha_radius_nroots(1);
D=k_D_alpha_radius_nroots(2);
alpha=k_D_alpha_radius_nroots(3);
radius=k_D_alpha_radius_nroots(4);
nroots=k_D_alpha_radius_nroots(5);
kaD=k*radius/D;
%determine the roots:
fun = @(x)(3*kaD.*x-alpha*x.^3)./(3*kaD+(kaD-1)*alpha*x.^2)-tan(x);
%% non zero positive roots up to 100 are searched:
xs = linspace(0.1,100,50001)';
ys = fun(xs);
scinter = find(diff(sign(ys)));
ninter = numel(scinter);
xroots = NaN(1,ninter);
%options = optimset('FunValCheck','on');
for i = 1:ninter
    xroots(i) = fzero(fun,xs(scinter(i) + [0 1]));%,options);
end
%% NB: for some reason fzero does not take into account the definition
%% domain of the function and returns zeros for values at which it is not defined.
%% The easy way to remove the "false positive" was to keep every second
%% root.
%for i = 1:1: nroots

    
for i = 1:1:size(xroots,2)/2
    q(i,1) = xroots(1,2*i);
end
end