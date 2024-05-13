function [pdf, logpdf] = halftpdf(x, nu, A)
% return a function proportional to the PDF of half student-t
% nu: degree of freedom

pdf=(1+(x./A).^2./nu).^(-(nu+1)./2);
logpdf=-(nu+1)./2*log1p((x./A).^2./nu);

end