function [g1] = get_g1(g, dim, TR)
% Return g1 unstandardized, computed by trapazoid rule or simpson's rule
% input:
%   g: generator, defined on [0, TR]
% output:
%   g1: spline funtion of g1

u=linspace(0, TR, 100)';
n=length(u);
g1_grid=zeros(n, 1);

% integrand for g1's evalution
const=(pi^((dim-1)/2) / gamma((dim-1)/2));
% Beware the integrand's last term is not defined at zero when dim <3
f_int=@(u, s) const*g(u+s).*s.^((dim-3)/2); 

% compute the g1(u) via trapzoid rule or Simpson's Rule
for j=1:(n-1)
    if dim==2
        X=linspace(10e-6, TR-u(j), 1000);
    else
        X=linspace(0, TR-u(j), 1000);
    end
    Y=f_int(u(j), X);
    Y(1)=Y(2);

    %g1_grid(j)=trapz(X, Y);
    g1_grid(j)=simps(X, Y);
end
g1_grid(n)=0;

% obtain smooth spline function form of g1
g1=fit(u, g1_grid, 'smoothingspline');


end