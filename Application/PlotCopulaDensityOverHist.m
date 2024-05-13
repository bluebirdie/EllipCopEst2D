%clearvars
close all


[u1,u2] = meshgrid(0.01:0.01:0.99);
[m, n]=size(u1);
UV=[u1(:), u2(:)];
Z_s_chain=zeros(m*n, nloop-nwarmup);

for i=1:(nloop-nwarmup)
    Weights=weights_chain(:, i);
    
    [g_s, f1, ~, Q1]= getFuncs(dim, Weights, Bdensities, xr, a, b);
    Z_s=EllipCopulaPDF(UV, Omega, g_s, f1, Q1);
    Z_s_chain(:, i)=Z_s;

    fprintf("loop: %d\n", i)
end
Z_s_mean=mean(Z_s_chain, 2);
Z_s=reshape(Z_s_mean, m, n);



figure
surfc(u1, u2, Z_s, 'FaceAlpha', 0.7)
%zlim([0,20])
hold on
%histogram2(U(:, 1), U(:,2), 'Normalization','pdf', 'FaceAlpha', 0.7, 'FaceColor','#EDB120',  'BinWidth',[0.15, 0.15])
histogram2(U(:, 1), U(:,2), 'Normalization','pdf', 'FaceAlpha', 0.7, 'FaceColor','#EDB120')


figure
histogram2(U(:, 1), U(:,2), 'Normalization','pdf', 'FaceColor','#EDB120')
zlim([0, 6])


ntt=1000;
tt=linspace(0, 1, ntt)';

g_chain=zeros(length(tt), nloop-nwarmup);
TR=((xr+a)^(dim/2)-a^(dim/2))^(2/dim);
for i=1:(nloop-nwarmup)
    Weights=weights_chain(:, i);
    g_est=get_g_Bspline_Lieb(dim, Weights, Bdensities, a);
    g_sest = g_standardize_Bdensity(g_est, TR, dim, b);
    g_chain(:, i)=g_sest(tt);

    fprintf("loop: %d\n", i)
end
g=@(t) g_Gaussian(t);
g_s=g_standardize_full(g, dim, b);
figure
plot(tt, mean(g_chain, 2))
hold on
plot(tt, g_s(tt), 'LineStyle','--')
legend("Estimated generator", "Generator of Gaussian copula")


figure
histogram(SnB)

