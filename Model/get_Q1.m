function [Q1] = get_Q1(cdf1, TR)
% return unstandardized Q1

ql=-sqrt(TR);
qr=sqrt(TR);

x=linspace(ql, qr, 500)';
y=cdf1(x);

y(1)=0;
y(end)=1;

[y_unique, ia] = unique(y);
x_unique=x(ia);

Q1=@(u) interp1(y_unique, x_unique, u);

end