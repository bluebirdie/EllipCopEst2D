function [Bdensities, T] = get_Bdensity(xl, xr, nbasis)
% create quadratic B-spline densities, modified from the Appendix of (Staudenmayer et al., 2008)
% Bdensities: cell array of basis densities
% T: knot sequence
% nbasis : number of basis densities, including 2 extras with different
% shapes at left side

% xl: left bound of knots
% xr: right bound of knots

% K: number of regular shaped basis densities in [xl, xr]
% the total number of basis functions is K+2, including 2 extras with different
% shapes at the left end

% each quadratic basis span over 4 knots
% so number of knots in [xl, xr] is K+3


% number of internal basis functions
K=nbasis-2;
% knots are equally spaced with delta
delta=(xr-xl)/(K+2);

% define knot sequence
T=xl: delta: xr;

% define basis densities, totally K+2 basis densities

B1toK=cell(K, 1);
for j=1:K
    B1toK{j}=@(x) (T(j)<=x & x< T(j+1)).*((x-T(j))./delta).^2./2./delta+...
        (T(j+1)<=x & x< T(j+2)).*(-((x-T(j+1))./delta).^2+(x-T(j+1))./delta+1/2)./delta+...
        (T(j+2)<=x & x< T(j+3)).*(1-(x-T(j+2))./delta).^2./2./delta;
end

Bhead=cell(2, 1);


% define B_(-1), span from t_1 to t_2 
Bhead{1}=@(x) (T(1)<=x & x< T(2)).*(1-(x-T(1))./delta).^2./2.*6./delta;

% B0, span from t_1 to t_3
Bhead{2}=@(x) (T(1)<=x & x< T(2)).*(-((x-T(1))/delta).^2+(x-T(1))./delta+1/2)*6/(5*delta)+...
        (T(2)<=x & x< T(3)).*(1-(x-T(2))./delta).^2/2*6/(5*delta);


Bdensities=cat(1, Bhead, B1toK);

end