function [g_s, B] = g_standardize_full(g, dim, b)
% return standardized copula generator, which is unique.
% base on Algorithm one in Derumigny and Fermanian 2022
% g take column input, return column output

I1=integral(@(s) s.^(dim/2-1).*g(s), 0, Inf);
I2=integral(@(s) s.^(dim/2-3/2).*g(s), 0, Inf);

% f_int1=@(s) g(s).*s.^(dim/2-1);
% f_int2=@(s) g(s).*s.^((dim-3)/2);
% 
% if dim==2
%     X=linspace(10e-6, TR, 1000);
% else
%     X=linspace(0, TR, 1000);
% end
% 
% Y1=f_int1(X);
% I1=simps(X, Y1);
% 
% Y2=f_int2(X);
% I2=simps(X, Y2);

SD=2*pi^(dim/2)/gamma(dim/2);
SDminus=2*pi^((dim-1)/2)/gamma((dim-1)/2);

B=(b*SD*I1/(SDminus*I2))^2;
A=2*B^(dim/2)/(SD*I1);

g_s=@(t) A*g(B*t);


end