function [sample] = simulation_MH(ndat, densityFUN, Delta)
% ndat: sample size
% densityFUN: target density
% Delta: tuning parameter
niter = 100;
sample = ones(ndat, 1);

for i = 1:ndat
    oldProb = densityFUN(sample(i));
    for iiter = 1:niter

        newPoint=unifrnd(sample(i)-Delta, sample(i)+Delta);
        newPoint=abs(newPoint);

        newProb = densityFUN(newPoint);

        if isfinite(newProb) && newProb > 0 && rand <= newProb/oldProb
            sample(i) = newPoint;
            oldProb = newProb;
        end
    end
end


end