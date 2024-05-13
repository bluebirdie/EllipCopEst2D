function [gt] = g_StudentT(t, dim, nu)

gt=(1+t./nu).^(-(nu+dim)/2);

end