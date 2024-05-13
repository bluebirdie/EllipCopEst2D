function X = sampleDist(f, M, N, b)
% SAMPLEDIST  Sample from an arbitrary distribution
%   f : target PDF
%   b : range, [min, max].
%   M : the threshold value for the proposal distribution, such that
%   f(x) < M for all x in b. ???
%   N: sample size
%
%   Examples:
%   Sample from a step function over [0,1]:
%   X = sampleDist(@(x)1.3*(x>=0&x<0.7)+0.3*(x>=0.7&x<=1),...
%                    1.3,1e6,[0,1]);
%   Sample from a normal distribution over [-5,5]:
%   X = sampleDist(@(x) 1/sqrt(2*pi) *exp(-x.^2/2),...
%                    1/sqrt(2*pi),1e6,[-5,5]);
%

% Dmitry Savransky (dsavrans@princeton.edu)
% May 11, 2010

n = 0;
X = zeros(N,1);
counter = 0; %???

while n < N && counter < 1000
    x = b(1) + rand(2*N,1)*diff(b); % uniform proposal, propose 2N candidates
    uM = M*rand(2*N,1);
    x = x(uM < f(x));
    if isempty(x)
        error('No Points Sampled!')
    end
    X(n+1:min([n+length(x),N])) = x(1:min([length(x),N - n]));
    n = n + length(x);
    counter = counter+1;
end

