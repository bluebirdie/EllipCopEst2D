function [gt] = g_sim4(t)

gt=exp(-t)+exp(-t/3).*cos(t).^2;

end