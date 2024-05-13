%% Plot trace

% Acceptance Rates
figure
plotTrace(e_chain)
title(strcat("Weights accepted by ", num2str(round(mean(e_chain), 2))))



% Weights
figure
for j=1:3
    subplot(1, 3, j)
    plotTrace(weights_chain(j, :))
    title(strcat("$w$", num2str(j)))
end
sgtitle('Trace Plot for Weights')

figure
for j=6:10
    subplot(1, 5, j-5)
    plotTrace(weights_chain(j, :))
    title(strcat("$w$", num2str(j)))
end
sgtitle('Trace Plot for Weights $6$ to $10$')

figure
plotTrace(tau_chain)
title("$\tau$")