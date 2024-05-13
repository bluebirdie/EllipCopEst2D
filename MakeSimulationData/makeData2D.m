clearvars
close all


ndat=3000;

dim=2;

b=1;


g_sim=@(t) g_sim1(t); 

dataname="Sim12D";

g_s=g_standardize_full(g_sim, dim, b);


% get PDF of Y= Psi(R^2)=-a+(a^{d/2}+R^{d})^{2/d}, where a=1
a=1;
fy = get_fy_from_g(dim, g_s, a);

xl=0;
xr=2;

yy=xl:0.01: xr;
figure
subplot(1, 2, 1)
plot(yy, fy(yy), 'Color', 'blue')
title("Standardized $f_{R^d}$")
subplot(1, 2, 2)
plot(0:0.01:1, g_s(0:0.01:1), 'Color', 'red')
title("$g_s$")

% get the true marginal pdf
f1=@(x) sd_EllipCop(dim-1)*integral(@(y) g_s(x.^2+y).*(y.^(dim/2-3/2)), 0, Inf);
f1vec=@(x) arrayfun(f1, x);
% get marginal CDF
cdf1=@(x) integral(f1vec, -Inf, x);
cdf1vec= @(x) arrayfun(cdf1, x);
% get sample of F1 on grid
F1x=linspace(-3, 3, 100)';
F1y=cdf1vec(F1x);
% get a spline function of cdf1, for more efficient computation later
cdf1_sp=fit(F1x, F1y, 'linearinterp');



M=max(fy(yy))+1;

%%% Make copula data
Omega=[1 0.2;
    0.2 1];
% Compute Cholesky factor A from Omega
A=chol(Omega, 'lower');

% the number of Replicates
nRep=1;

for iRep=1:nRep
    % Sample Y by rejection sampling
    Y = sampleDist(@(x) fy(x), M, ndat, [xl, xr]);
    % Generate Usp from sphere
    Usp=unispherernd(dim, ndat);
    % compute R
    a=1;
    R=((Y+a).^(dim/2)-a^(dim/2)).^(1/dim);
    % Compute Z=RAU
    Z=R.*(A*Usp')';

    % get copula Data: U=F1(X)
    U=reshape(cdf1_sp(Z), size(Z));
    filename=strcat(dataname, num2str(iRep));
    save(filename, "U", "dim", "ndat", "Omega", "g_s", "b", "xr");

end