function[q] = Mizusaki_D_roots(alpha_nroots)
%Determine the roots of tan(x)=3*x/(3+alpha*x^2)
sym x;
alpha=alpha_nroots(1);
nroots=alpha_nroots(2);
%determine the roots:
fun = @(x) 3*x./(3+alpha*x.^2)-tan(x);
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
for i = 1:1:nroots
    q(i,1) = xroots(1,2*i);
end
end