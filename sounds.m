function [] = sounds(num, volume)   % Put in an amp/loudness parameter

if nargin < 2, volume = 1; end

switch num
    case -1
        fs = 20500;
        duration = 0.2;    % Measured in seconds
        freq = 100;
        values = 0:1/fs:duration;
        a = volume*sin(2*pi* freq*values);
        sound(a)
    case 0
        w = 100;
        t = linspace(0,2*pi,1800);
        y = volume*cos(w*t);
        sound(y);
    case 1     % Beep for Correct
        w = 100;
        t = linspace(0,2*pi,900);
        y = volume*cos(w*t);
        sound(y);
    case 2 % separate sound for 'no response': three rapid pulses of the 'incorrect' pitch
        w = 100;
        t = linspace(0,2*pi,2400);
        y = volume*cos(w*t);
        y(700:900) = 0;
        y(1500:1700) = 0;
        sound(y);
        
    otherwise
        error('invalid sounds number');
end
end