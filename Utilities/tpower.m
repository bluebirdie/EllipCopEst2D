function [tp] = tpower(x, knot, p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tp = (x - knot).^p.*(x > knot);

end