function [] = plotCred(X, L, U)
% plot credible interval
% X: x axis
% L: lower bound
% U: upper bound


if size(X, 1)==1
    % for row vectors
    patch_x=[X, fliplr(X)];
    patch_y=[L,  fliplr(U)];
    patch(patch_x, patch_y, [0.9290 0.6940 0.1250], ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none')
else
    % for column vectors
    patch_x=[X; flipud(X)];
    patch_y=[L; flipud(U)];
    patch(patch_x, patch_y, [0.9290 0.6940 0.1250], ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none')
end


end