function [Weights] = Etas_to_Weights(Etas)
% unconstrained variable Etas to constrained mixing proportions
% Etas: log(w_j/w_K), K-1 elements, because Eta_K=0
% Weights: K elements, 0 to 1 and sum to 1

% Transform back to constrained, Weights = exp(Etas)/sum(exp(Etas))

EtasFull=[Etas; 0];
Weights=softmax(EtasFull);

end