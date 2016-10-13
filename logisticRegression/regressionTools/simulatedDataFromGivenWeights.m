function[simulated_choice]=simulatedDataFromGivenWeights(fixed_weights,obtained_clicks)

simulated_choice = fixed_weights'*obtained_clicks;
for i=1:length(choice)
    if simulated_choice(i)>=0
        simulated_choice(i)=1;
    else
        simulated_choice(i)=0;
    end
end
end
