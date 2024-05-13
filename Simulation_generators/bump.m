function [out] = bump(x)
% a smooth function supported on [1, 1+pi]

out=(x>=1).*(x<=(1+pi)).*(x-1).*(1+pi-x).*sin(x-1);

end