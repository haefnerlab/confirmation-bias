function [order_of_orientations, correct_orientation_answer] = makeOrientations(desired_orientations, probabilities, number_of_images)

% makeOrientations is a function to generate a vector of orientation angles for the gabor patches

order_of_orientations = zeros(1,number_of_images);
    % an array to say what order of orientations the images should be presented in

different_probabilities = zeros(1,length(probabilities));
    % an array to say what the dividing probability values to have
    % Ex. If probabilities is [.1,.4,.3,.2], then the for loop below will
    % produce [.1,.5,.8,1]. This allows us to see which value rand will
    % fall between and determine which angle to have for the current image
    
    % Note that different_probabilities always has the same size as
    % desired_orientations and probabilities

for p = 1:length(different_probabilities)
    current_probability = 0;
    for k = 1:p
        current_probability = current_probability + probabilities(k);
    end
    % This inner for loop adds up all of the probability values in the
    % array, probabilities, which has an index smaller than p as p iterates
    % through the array, different_probabilities
    
    different_probabilities(p) = current_probability;
end


for i = 1:number_of_images
    random_number = rand;
    for j = 1:length(different_probabilities)
        if random_number < different_probabilities(j)
            order_of_orientations(i) = desired_orientations(j);
            break;
        end
    end
    % This inner for loop iterates through different_probabilities to check
    % if random_number is smaller than the current probability value being
    % examined. Ex. If different_probabilities = [.1,.5,.8,1] and rand =
    % .05, then since .05 < .1, whatever orientation angle is in the first
    % index of desired_orientations is choosen. Otherwise if rand = .92,
    % then the for loop would iterate to the last index of
    % different_probabilities before seeing .92 < 1 and choosing the
    % orientation in the last index of desired_orientations.
    
end

correct_orientation_answer = mode(order_of_orientations);
        % Subjects are expected to be able to chose the most frequent orientation
        
        % If the most frequent orientations have the exact same frequency, mode will output
        % the smaller orientation value
end