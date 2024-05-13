function [Usp]=unispherernd(dim, ndat)
% dim=2 or 3

if dim==3
    u1=rand(ndat,1);
    u2=rand(ndat,1);

    lambda=acos(2*u1-1)-0.5*pi;
    phi=2*pi*u2;

    x=cos(lambda).*cos(phi);
    y=cos(lambda).*sin(phi);
    z=sin(lambda);

    Usp=[x, y, z];
else
    theta=2*pi.*rand(ndat, 1);
    Usp=[cos(theta), sin(theta)]; % theta is radian
end

end