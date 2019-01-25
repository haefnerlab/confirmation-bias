function c = betacolor(beta, redmax, blumax)
%MODEL.BETACOLOR helper function to for colormap of model slope (beta). Inputs are beta, the minimum
%expected beta (redmax < 0) and maximum expected beta (blumax > 0). 
blu = [0 0 1];
red = [1 0 0];
zero = [.6 0 .6];

if (redmax == 0 && beta < 0) || (blumax == 0 && beta > 0)
    % Special case when blumax or redmax is set to zero, don't allow any blue (or red respectively)
    % colors at all.
    c = zero;
elseif beta < 0
    betafrac = min(1, beta / redmax);
    c = red * betafrac + zero * (1 - betafrac);
else
    betafrac = min(1, beta / blumax);
    c = blu * betafrac + zero * (1 - betafrac);
end