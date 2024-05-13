clearvars
close all


load("Sim1.mat");


% the number of iterations in MCMC loop
nloop=10000;
nwarmup=round(0.5*nloop); % number of iterations used as "burn-in" stage


% Tuning parameters
Eta_tuning=0.1;
tau_trans_tuning=0.5;

% the number of basis densities
nbasis=30;
xr=1.5; % basis densities will be constructed in interval [0, xr]

[Bdensities, Knots] = get_Bdensity(0, xr, nbasis);

a=1; % a positive number in Liebscher's transormation function
b=1; % a positive number for standardizing the generator




% Initialize parameters with starting values
% etas are nu's in the paper, which are uncontrained parameters transformed
% from weights
Etas=zeros(nbasis-1, 1);
% obtain weights from etas
Weights=Etas_to_Weights(Etas);
% tau_trans is the tau_star in the paper
tau_trans=normrnd(0, 0.1);

% hold values after burn-in stage, each column is one sample from an
% iteration
weights_chain=zeros(nbasis, nloop-nwarmup);
tau_chain=zeros(1, nloop-nwarmup);
e_chain=zeros(1, nloop-nwarmup);

for i=1:nloop

    if mod(i, 1)==0
        fprintf("loop: %d\n", i)
    end

    % Draw Etas
    Etas1=unifrnd(Etas-Eta_tuning, Etas+Eta_tuning);
    % Draw tau_trans
    tau_trans1=unifrnd(tau_trans-tau_trans_tuning, tau_trans+tau_trans_tuning);

    % Metropolis steps
    A=logPostCopulaBspline_Lieb(U, dim, Omega, Bdensities, Etas1, tau_trans1, xr, a, b);
    B=logPostCopulaBspline_Lieb(U, dim, Omega, ...
        Bdensities, Etas, tau_trans, xr, a, b);

    logRate=A-B;
    current_loglike=B;
    Rate=exp(logRate);
    e=min(1,Rate);
    u=rand;

    if u<=e
        Etas=Etas1; %update eta
        tau_trans=tau_trans1; %update tau_trans
        Weights=Etas_to_Weights(Etas);
    end

    % record value after burn in
    % each COLUMN is a sample
    if i>nwarmup
        weights_chain(:, i-nwarmup)=Weights;
        tau_chain(i-nwarmup)=tau_trans;
        e_chain(i-nwarmup)=e;
    end

end



%% Plot 
setupPlot;


ntt=1000;


tt=linspace(0, 1, ntt)';


TR=((xr+a)^(dim/2)-a^(dim/2))^(2/dim);
g_chain=zeros(length(tt), nloop-nwarmup);

for i=1:(nloop-nwarmup)
    Weights=weights_chain(:, i);
    g_est=get_g(dim, Weights, Bdensities, a);
    g_sest = g_standardize(g_est, TR, dim, b);
    g_chain(:, i)=g_sest(tt);
end


g_mean=mean(g_chain, 2);


figure
plot(tt, g_mean,'Color', 'Blue', 'LineStyle','-')
hold on
plot(tt, g_s(tt), "Color", 'black', 'LineStyle','--')
xlabel("$t$")
ylabel("$\hat g(t)$")
legend("$\hat g$", "True $g$")
ax = gca; % get current axes
ax.TickLabelInterpreter = 'latex';
title("Estimated generator and the true one")

savefig("fitted_generator")



%% Save results

save("Simulation1Result.mat");




