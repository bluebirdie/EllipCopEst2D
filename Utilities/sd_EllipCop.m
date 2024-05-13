function [output] = sd_EllipCop(dim)
% return normalizing constant of PDF of R^2

output=(pi^(dim/2))/gamma(dim/2);
end