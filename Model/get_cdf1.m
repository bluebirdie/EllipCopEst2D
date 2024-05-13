function [cdf1] = get_cdf1(f1_sp, TR)
% f1_sp: unstandardized f1 as spline functiom form
% cdf1:unstandardized F1

xll=-sqrt(TR);
cdf1=@(x) min(1, integrate(f1_sp, x, xll));

end