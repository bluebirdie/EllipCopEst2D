function [f1] = get_f1(g1, TR)
% return f1(unstandardized)
% g1: unstandardized g1

f1x=linspace(-sqrt(TR), sqrt(TR), 100)';
f1y=g1(f1x.^2);

% unnormalized f1
f10=fit(f1x, f1y, 'smoothingspline');
% calculate normalizing constant of f10
n_constant=integrate(f10, sqrt(TR), -sqrt(TR));
% obtain coefficients of cfit
coeffvals= coeffvalues(f10);
% update coefficients of cfit to reinforce normalization
coeffvals.coefs=coeffvals.coefs./n_constant;

f1=cfit(fittype('smoothingspline'), coeffvals);


end

