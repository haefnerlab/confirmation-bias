function c = betacolor(beta, redmax, blumax)
blu = [0 0 1];
red = [1 0 0];
zero = [.6 0 .6];

if beta < 0
    betafrac = min(1, beta / redmax);
    c = red * betafrac + zero * (1 - betafrac);
else
    betafrac = min(1, beta / blumax);
    c = blu * betafrac + zero * (1 - betafrac);
end