function[]=TestAnalysisAuditory()

%% Make clicks

number_of_trials = 1000;
leftVolume=0.1;
rightVolume = 0.1;
number_of_left_clicks=zeros(number_of_trials);
number_of_right_clicks = zeros(number_of_trials);
bins=120;
clicks=zeros(number_of_trials,2,bins);
left_click_rate=20;
right_click_rate=4;
for i=1:number_of_trials
    [number_of_left_clicks(i), number_of_right_clicks(i), clicks(i,:,:)] = makePoissonClicks(bins, left_click_rate, right_click_rate, leftVolume, rightVolume);
end
order_of_clicks = [squeeze(Test_Data.clicks(:,1,:)) squeeze(Test_Data.clicks(:,2,:))];
order_of_clicks = reshape(order_of_clicks, number_of_trials, []);
simulated_choice = simulatedDataFromGivenWeights(fixed_weights,order_of_clicks)

end