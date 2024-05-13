function [] = plotTrace(X)
% Create Trace plot of X
% X: row or colum vector

n=length(X);
scatter(1:n, X, 10, 'MarkerFaceColor', "#0072BD", 'MarkerEdgeColor','none');

end