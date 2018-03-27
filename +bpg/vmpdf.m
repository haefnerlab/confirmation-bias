function p = vmpdf(x, theta, kappa)
%VMPDF Von-Mises circular pdf. Theta is mean, kappa controls width. Inputs
%are in radians.
p = exp(kappa * cos(x - theta)) / besseli(0, kappa);
end