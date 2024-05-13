function [f] = BdensityMix(x, Weights, Bdensities)
% return a Bdensity mixture 

f=0;
for j=1:length(Weights)
   f=f+ Weights(j)*Bdensities{j}(x);
end

end