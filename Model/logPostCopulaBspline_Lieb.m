function [logPost, g, f1, cdf1, Q1, Rsq] = logPostCopulaBspline_Lieb(U, dim, Omega, Bdensities, Etas, tau_trans, xr, a, b)
% Return log posterior PDF (unormalized) of parameters

% U: Copula data ndat by dim

% a: positive constant for Liebcher transformation
% b: positive constant for generator standardization

ndat=size(U,1);

% Obtain weights from Etas
Weights=Etas_to_Weights(Etas);
nbasis=length(Weights);
tau=exp(tau_trans);

% Obtain related functions
[g, f1, cdf1, Q1]= getFuncs(dim, Weights, Bdensities, xr, a, b);


% Obtain pseudo data, which has elliptical distribution
Z=reshape(Q1(U(:)), size(U));

% Obtain R^2 in Stochastic form
Rsq=zeros(ndat, 1); % ndat by 1 vector
for i=1:ndat
    Rsq(i)=Z(i, :)/Omega*Z(i, :)'; 
end

% when values are too small <10e-6, impute them with 10e-6 to avoid
% round-off error
logNum=reallog(max(g(Rsq(:)),10e-6));
logDen=reallog(max(f1(Z(:)),10e-6)); 

logLike=sum(logNum, 'all')-sum(logDen, 'all');

% Obtain log prior of etas

% second order difference matrix
D = diff(eye(nbasis-1), 2);
P = D'*D;
c=2;
P(1, 1)=P(1, 1)+1/c;
P(2, 2)=P(2, 2)+1/c;

logPEtas=(1-nbasis)*log(tau)-1/(2*tau^2)*(Etas'*P*Etas);

nu_tau=5;
A_tau=100;

% log prior of tau
[~, logPtau] = halftpdf(tau, nu_tau, A_tau);
% Jacobian for transform from tau to tau_trans
logJtau=log(tau);

logPost=logLike+logPtau+logJtau+logPEtas;

end