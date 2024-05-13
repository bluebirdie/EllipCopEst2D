function [g] = get_g(dim, Weights, Bdensities, a)
% Return the generator (unstandardized) based on a B-spline Density Mixture
% defined in [0, TR=((XR+a)^d/2-a^d/2)^2/d]
% Weights: mixing proportions
% a: positive constant for Liebcher transformation

const=gamma(dim/2)/pi^(dim/2);

g=@(t) const*(a^(dim/2)+t.^(dim/2)).^(2/dim-1).*...
    BdensityMix(max(-a+(a^(dim/2)+t.^(dim/2)).^(2/dim), 0), ...
    Weights, Bdensities);

end