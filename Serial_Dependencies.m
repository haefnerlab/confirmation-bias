function [prob_correct_left, prob_correct_right, prob_wrong_left, prob_wrong_right, correct_left_trials, correct_right_trials, wrong_left_trials, wrong_right_trials] = Serial_Dependencies(Data)

% Example Input:
%[LcL, LcR, RcL, RcR, LaL, LaR, RaL, RaR] = Serial_Dependencies(Test_Data, -pi/4, pi/4)

trials = length(Data.choice);


prob_correct_left = 0;    % Subject chose left correctly on the prev trial
prob_correct_right = 0;   % Subject chose right correctly on the prev trial
prob_wrong_left = 0;      % Subject chose left incorrectly on the prev trial
prob_wrong_right = 0;     % Subject chose right incorrectly on the prev trial

correct_left_trials = 0;
correct_right_trials = 0;
wrong_left_trials = 0;
wrong_right_trials = 0;


for t = 1:trials
    if t ~= 1
        if Data.choice(t-1) == 1 && Data.accuracy(t-1) == 1
            % Previous trial choice = left and was accurate
            correct_left_trials = correct_left_trials + 1;
            if Data.choice(t) == 1
                % Does the subject chose left this time?
                prob_correct_left = prob_correct_left + 1;
            end
            
        elseif Data.choice(t-1) == 0 && Data.accuracy(t-1) == 1
            % Previous trial choice = right and was accurate
            correct_right_trials = correct_right_trials + 1;
            if Data.choice(t) == 1
                % Does the subject chose left this time?
                prob_correct_right = prob_correct_right + 1;
            end
            
        elseif Data.choice(t-1) == 1 && Data.accuracy(t-1) == 0
            % Previous trial choice = left and was inaccurate
            wrong_left_trials = wrong_left_trials + 1;
            if Data.choice(t) == 1
                % Does the subject chose left this time?
                prob_wrong_left = prob_wrong_left + 1;
            end
            
        elseif Data.choice(t-1) == 0 && Data.accuracy(t-1) == 0
            % Previous trial choice = right and was inaccurate
            wrong_right_trials = wrong_right_trials + 1;
            if Data.choice(t) == 1
                % Does the subject chose left this time?
                prob_wrong_right = prob_wrong_right + 1;
            end
            
        end
    end
end


end