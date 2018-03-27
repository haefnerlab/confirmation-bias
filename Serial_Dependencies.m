function [leftAfterCorrectLeft, leftAfterCorrectRight, leftAfterWrongLeft, leftAfterWrongRight] = Serial_Dependencies(Data)
choseLeft = Data.choice(1:end-1) == +1;
choseRight = ~choseLeft;

wasCorrect = Data.accuracy(1:end-1);
wasWrong = ~wasCorrect;

correctLeft = wasCorrect & choseLeft;
correctRight = wasCorrect & choseRight;
wrongLeft = wasWrong & choseLeft;
wrongRight = wasWrong & choseRight;

followedLeft = Data.choice(2:end) == +1;

leftAfterCorrectLeft  = sum(correctLeft  & followedLeft) / sum(correctLeft);
leftAfterCorrectRight = sum(correctRight & followedLeft) / sum(correctRight);
leftAfterWrongLeft    = sum(wrongLeft    & followedLeft) / sum(wrongLeft);
leftAfterWrongRight   = sum(wrongRight   & followedLeft) / sum(wrongRight);
end