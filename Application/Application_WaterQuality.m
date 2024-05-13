clearvars
close all


load('CopulaData.mat')
N=size(U, 1);
dim=size(U, 2);



nloop=10000;
nwarmup=5000;

% Tuning parameters
Eta_tuning=0.1;
tau_trans_tuning=0.5;

nbasis=20;
a=1;
b=1;
xr=2;
TR=((xr+a)^(dim/2)-a^(dim/2))^(2/dim);
[Bdensities, Knots] = get_Bdensity(0, xr, nbasis);

figure
xx=0:0.01:xr;
hold on
for iplot=1:nbasis
    plot(xx, Bdensities{iplot}(xx))
end


%% Initialize parameters
% Weights trans
Etas=zeros(nbasis-1, 1);
%Etas=normrnd(0, 0.1, nbasis-1, 1);
Weights=Etas_to_Weights(Etas);
tau_trans=normrnd(0, 0.1);

% hold values after burn in, each column is one sample
weights_chain=zeros(nbasis, nloop-nwarmup);
tau_chain=zeros(1, nloop-nwarmup);
e_chain=zeros(1, nloop-nwarmup);
SnB=zeros(1, nloop-nwarmup);
for i=1:nloop

    if mod(i, 1)==0
        fprintf("loop: %d\n", i)
    end

    % Draw Etas
    Etas1=unifrnd(Etas-Eta_tuning, Etas+Eta_tuning);

    % Draw taus
    tau_trans1=unifrnd(tau_trans-tau_trans_tuning, tau_trans+tau_trans_tuning);

    A=logPostCopulaBspline_Lieb(U, dim, Omega, Bdensities, Etas1, tau_trans1, xr, a, b);
    B=logPostCopulaBspline_Lieb(U, dim, Omega, Bdensities, Etas, tau_trans, xr, a, b);
    % disp(B)
    logRate=A-B;
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

