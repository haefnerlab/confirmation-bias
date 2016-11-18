function GaborData = newGaborData(varargin)

    function value = get_arg(name, default)
        % Helper function to get named arguments with a default
        idx = strcmpi(name, varargin);
        if any(idx)
            val_idx = find(idx)+1;
            value = varargin{val_idx};
            varargin(find(idx):val_idx) = [];
        else
            value = default;
        end
    end

%% User-settable params
GaborData.trials_per_block = get_arg('trials_per_block', 100);
GaborData.blocks = get_arg('blocks', 4);
GaborData.stair_fn = get_arg('stair_fn', @GaborStaircase.contrast);

total_trials = GaborData.trials_per_block * GaborData.blocks;

% Initial values of staircase-able parameters
GaborData.contrast = zeros(1, total_trials);
GaborData.contrast(1) = get_arg('contrast', 32);
GaborData.ratio = zeros(1, total_trials);
GaborData.ratio(1) = get_arg('ratio', 0.8);
GaborData.pixel_noise = zeros(1, total_trials);
GaborData.pixel_noise(1) = get_arg('pixel_noise', 1); % standard deviation of pixel noise
GaborData.step_size = zeros(1, total_trials);

% Staircase bounds and step size, with defaults set depending on stair_fn
GaborData.model_observer = get_arg('model_observer', '');
if isequal(GaborData.stair_fn, @GaborStaircase.contrast)
    GaborData.stair_bounds = get_arg('stair_bounds', [0 64]);
    GaborData.step_size(1) = get_arg('step_size', 2); % multiplicative (in the "easier" direction)
elseif isequal(GaborData.stair_fn, @GaborStaircase.ratio)
    GaborData.stair_bounds = get_arg('stair_bounds', [0.5 1.0]);
    GaborData.step_size(1) = get_arg('step_size', 1); % additive (in the "easier" direction)
elseif isequal(GaborData.stair_fn, @GaborStaircase.pixel_noise)
    GaborData.stair_bounds = get_arg('stair_bounds', [0 32]);
    GaborData.step_size(1) = get_arg('step_size', -4); % additive (in the "easier" direction)
end

% Other misc. user-definable parameters
GaborData.number_of_images = get_arg('number_of_images', 10);
GaborData.stimulus_fps = get_arg('stimulus_fps', 12);	% frame rate of stimuli
GaborData.screen_resolution = get_arg('screen_resolution', 25);          % how many pixels correspond to a single datapoint of a gabor
GaborData.image_length_x = get_arg('image_length_x', 5);  % Size of the image along x-axis
GaborData.image_length_y = get_arg('image_length_y', 5);

% Preallocate fields that will be populated with data by running the
% experiment.
GaborData.streak = zeros(1, total_trials);
GaborData.reversal_counter = zeros(1, total_trials);
GaborData.correct_answer = zeros(1, total_trials);
GaborData.ideal_answer = zeros(1, total_trials);
GaborData.reaction_time = zeros(1, total_trials);
GaborData.choice = zeros(1, total_trials);
GaborData.accuracy = zeros(1, total_trials);
GaborData.order_of_orientations = zeros(total_trials, GaborData.number_of_images);
GaborData.log_frame_odds = zeros(total_trials, GaborData.number_of_images);
GaborData.log_decision_odds = zeros(total_trials, GaborData.number_of_images);
GaborData.average_orientations = zeros(2, total_trials);

GaborData.current_trial = 0;

GaborData.eye_tracker_points = {};

GaborData.left_template = eye(GaborData.image_length_y, GaborData.image_length_x);
GaborData.right_template = rot90(GaborData.left_template);

if ~isempty(varargin)
    warning('Unkown arguments given to newGaborParams');
end

% Sanity checks for common "gotchas"
if ~isempty(GaborData.model_observer) && ~isempty(GaborData.stair_fn)
    warning('Model observer with a staircase?');
end

if isequal(GaborData.stair_fn, @GaborStaircase.ratio)
    if GaborData.step_size(1) < 0
        warning('Changing sign of ratio step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
end

if isequal(GaborData.stair_fn, @GaborStaircase.pixel_noise)
    if GaborData.step_size(1) > 0
        warning('Changing sign of pixel_noise step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
end

if isequal(GaborData.stair_fn, @GaborStaircase.contrast)
    if GaborData.step_size(1) < 0
        error('Contrast staircase is multiplicative; step size of %f doesn''t make sense', GaborData.step_size(1));
    elseif GaborData.step_size(1) < 1
        warning('Chaning contrast step_size < 1 to 1/step_size');
        GaborData.step_size(1) = 1 / GaborData.step_size(1);
    end
end

if GaborData.ratio(1) > 1 || GaborData.ratio(1) < 0
    error('Ratio should be between 0 and 1');
end

if ~isempty(GaborData.model_observer) && ~any(strcmpi(GaborData.model_observer), {'ideal'})
    warning('%s is not a known model observer', GaborData.model_observer);
end

disp(GaborData);

end