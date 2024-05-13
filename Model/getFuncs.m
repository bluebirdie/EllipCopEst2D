function [g_s, f1_s, cdf1_s, Q1_s] = getFuncs(dim, Weights, Bdensities, xr, a, b)
% Return standardized functions according to  parameters of Bdensity
% mixture

TR=((xr+a)^(dim/2)-a^(dim/2))^(2/dim);

% g (unstandardized), defined on [0, TR=((XR+a)^d/2-a^d/2)^2/d]
g=get_g(dim, Weights, Bdensities, a);
% obtain standardize generator, denoted as g_s
[g_s, A] = g_standardize(g, TR, dim, b);


% g1 (unstandardized), defined on [0, TR=xr^(2/dim)]
g1 = get_g1(g, dim, TR);


% f1 (unstandardized), defined on [-sqrt(TR), sqrt(TR)]
f1 = get_f1(g1, TR);
% f1(standardized), x must be column vector
f1_s=@(x) (sqrt(A)*f1(sqrt(A)*x)).*(x>=-sqrt(TR/A) & x<=sqrt(TR/A));



% F1 (unstandardized), defined on [-sqrt(TR), sqrt(TR)]
cdf1 = get_cdf1(f1, TR);
% F1 (standardized)
cdf1_s=@(x) ( x<=-sqrt(TR/A) )*0+...
    ( x>-sqrt(TR/A) & x<sqrt(TR/A) ).*cdf1(sqrt(A)*x)+...
    ( x>=sqrt(TR/A) )*1;


% marginal quantile Q1 (unstandardized), with range on [-sqrt(TR), sqrt(TR)]
Q1 = get_Q1(cdf1, TR); 
% Q1 (standardized), with range on [-sqrt(TR)/sqrt(A), sqrt(TR)/sqrt(A)]
Q1_s=@(p) Q1(p)/sqrt(A);



end