function [g_s, A] = g_standardize(g, TR, dim, b)
% return standardized copula generator g_s, which is unique.
% g is a normalized generator
% g take column input, return column output


if dim==2 %when dimension is 2, the integrand is not defined at X=0, see line 13
    X=linspace(10e-5, TR, 1000);
else
    X=linspace(0, TR, 1000);
end

Y=X.^((dim-3)/2).*g(X);
% Choose trazoidal rule or Simpson's rule
%I2=trapz(X, Y);
I2=simps(X, Y);

% standardize g
A=(b/I2/sd_EllipCop(dim-1))^2;
g_s=@(t) (t<=TR/A).*A^(dim/2).*g(A.*t);

end