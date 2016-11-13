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
GaborData.staircase = get_arg('staircase', 'contrast'); % contrast, ratio, or pixel_noise

total_trials = GaborData.trials_per_block * GaborData.blocks;

GaborData.contrast = zeros(1, total_trials);
GaborData.contrast(1) = get_arg('contrast', 64);
GaborData.ratio = zeros(1, total_trials);
GaborData.ratio(1) = get_arg('ratio', 0.85);
GaborData.pixel_noise = zeros(1, total_trials);
GaborData.pixel_noise(1) = get_arg('pixel_noise', 4); % standard deviation of pixel noise
GaborData.step_size = zeros(1, total_trials);

GaborData.model_observer = get_arg('model_observer', '');
if strcmpi(GaborData.staircase, 'contrast')
    GaborData.stair_bounds = get_arg('stair_bounds', [0 128]);
    GaborData.step_type = get_arg('step_type', 'add'); % add or multiply
    GaborData.step_size(1) = get_arg('step_size', -2); % sign should be in the "easier" direction
elseif strcmpi(GaborData.staircase, 'ratio')
    GaborData.stair_bounds = get_arg('stair_bounds', [0.5 1]);
    GaborData.step_type = get_arg('step_type', 'multiply'); % add or multiply
    GaborData.step_size(1) = get_arg('step_size', 1.1); % sign should be in the "easier" direction
end
GaborData.number_of_images = get_arg('number_of_images', 10);
GaborData.stimulus_fps = get_arg('stimulus_fps', 12);	% frame rate of stimuli
GaborData.screen_resolution = get_arg('screen_resolution', 25);          % how many pixels correspond to a single datapoint of a gabor
GaborData.image_length_x = get_arg('image_length_x', 5);  % Size of the image along x-axis
GaborData.image_length_y = get_arg('image_length_y', 5);

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
if ~isempty(GaborData.model_observer) && GaborData.step_size(1) ~= 0
    warning('Model observer with nonzero step size?');
end

if ~any(strcmpi(GaborData.step_type, {'add', 'multiply'}))
    error('Unknown step type, %s', GaborData.step_type);
end

if strcmpi(GaborData.staircase, 'contrast')
    if GaborData.step_size(1) > 0
        warning('Changing sign of contrast step_size from %d to %d', GaborData.step_size, -GaborData.step_size);
        GaborData.step_size = -GaborData.step_size;
    elseif strcmpi(GaborData.step_type, 'multiply')
        warning('Contrast step_type should be additive');
    end
end

if strcmpi(GaborData.staircase, 'ratio')
    if GaborData.step_size(1) < 0 || strcmpi(GaborData.step_type, 'add')
        warning('Ratio step_type should be multiplicative');
    elseif GaborData.step_size(1) < 1
        warning('Chaning step_size < 1 to 1/step_size');
        GaborData.step_size(1) = 1 / GaborData.step_size(1);
    end
end

if strcmpi(GaborData.staircase, 'pixel_noise') && GaborData.step_size(1) > 0
    warning('Changing sign of pixel_noise step_size from %d to %d', GaborData.step_size, -GaborData.step_size);
    GaborData.step_size = -GaborData.step_size;
end

if ~strcmpi(GaborData.staircase, 'contrast') && all(GaborData.stair_bounds == [0 128])
    warning('Need new stair_bounds for non-contrast staircase');
end

if GaborData.ratio(1) > 1 || GaborData.ratio(1) < 0
    error('Ratio should be between 0 and 1');
end

if ~isempty(GaborData.model_observer) && ~any(strcmpi(GaborData.model_observer), {'ideal'})
    warning('%s is not a known model observer', GaborData.model_observer);
end

disp(GaborData);

end