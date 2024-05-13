function [pdf] = EllipCopulaPDF(U, Omega, g, f1, Q1)
% Return a 2D copula density from mixParameters

% U: a ROW vector or ndat by dim matrix
% g: normalized generator
% f1: marginal pdf
% Q1: marginal quantile function

ndat=size(U,1);
Z=Q1(U); % a row or n by dim matrix

Rsq=zeros(ndat, 1); % ndat by 1 vector
for i=1:ndat
    Rsq(i)=Z(i, :)/Omega*Z(i, :)'; %make this more efficient!
end


pdf=det(Omega)^(-0.5)*g(Rsq)./prod([f1(Z(:, 1)), f1(Z(:, 2))], 2);


end