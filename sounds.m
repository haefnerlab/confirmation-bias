function [] = sounds(num)   % Put in an amp/loudness parameter

switch num
    case -1
        amp = 1.5;
        fs = 20500;
        duration = 0.2;    % Measured in seconds
        freq = 100;
        values = 0:1/fs:duration;
        a = amp*sin(2*pi* freq*values);
        sound(a)
    case 0
        w = 100;
        t = linspace(0,2*pi,1800);
        y = .2*cos(w*t);
        sound(y);
    case 1     % Beep for Correct
        w = 100;
        t = linspace(0,2*pi,900);
        y = .2*cos(w*t);
        sound(y);
   
    otherwise
        error('invalid sounds number');
    end
end