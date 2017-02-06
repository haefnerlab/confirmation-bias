function varargout = LoadOrRun(func, args, save_file, varargin)
%LOADORRUN load precomputed results from a file, or compute and save
%them to the file.
%
% ... = LOADORRUN(func, {arg1, arg2, ..}, save_file, ..) looks for results
% in save_file (a path to a .mat file). If the file does not exist,
% computes func(arg1, arg2, ..) and saves those results to the file. Return
% values are identical to whatever func returns.
%
% ... = LOADORRUN(..., '-recompute') forces results to be recomputed.
%
% ... = LOADORRUN(..., '-verbose') prints useful information.
%
% S = LOADORRUN(..., '-struct', 'fieldxyz') (only valid if the func returns
% a single struct, S). Use this form when func adds a single field
% 'fieldxyz' to a struct. Only the absence of the field S.fieldxyz will
% trigger the evaluation of func.

%% Parse extra flags and remove them from varargin
flag_recompute = strcmpi('-recompute', varargin);
varargin = varargin(~flag_recompute);
recompute = any(flag_recompute);

flag_verbose = strcmpi('-verbose', varargin);
varargin = varargin(~flag_verbose);
verbose = any(flag_verbose);

flag_struct = strcmpi('-struct', varargin);
struct_mode = any(flag_struct);
if struct_mode
    % 'right shift' to get logical array that indexes the next argument
    flag_field = [false flag_struct(1:end-1)];
    fieldname = varargin{flag_field};
    % clear both '-struct' and the fieldname from varargin
    varargin = varargin(~(flag_struct | flag_field));
    
    if nargout > 1
        error('Only 1 return value is allowed with the ''-struct'' flag');
    end
end

if ~isempty(varargin), warning('LoadOrRun got unexpected args.'); end

results = cell(1, nargout);

%% Ensure that parent directory of save_file exists.
save_dir = fileparts(save_file);
if ~isempty(save_dir) && ~exist(save_dir, 'dir')
    if verbose
        fprintf('Creating output directory: %s\n', save_dir);
    end
    mkdir(save_dir);
end

%% Determine whether a call to func is needed.
do_compute = ~exist(save_file, 'file') || recompute;

if ~do_compute && struct_mode
    % If file exists but the output doesn't contain 'fieldxyz' (in 'struct
    % mode'), then we still need to set do_compute=true
    contents = load(save_file);
    if ~isfield(contents.varargout{1}, fieldname)
        do_compute = true;
        if verbose
            fprintf('Computing struct field ''%s''\n', fieldname);
        end
    end
end

%% Call func or load precomputed results.
if do_compute
    % Call func(args) and capture as many return values as have been
    % requested by whoever called this function.
    if verbose
        fprintf('Calling %s with %d outputs...\t', func2str(func), nargout);
    end
    [results{:}] = func(args{:});
    
    % Sanity-check outputs when in struct mode
    if struct_mode
        if ~isstruct(results{1})
            error('func must return a struct when the ''-struct'' flag is used');
        elseif ~isfield(results{1}, fieldname)
            warning('In ''-struct'' mode, but func failed to add a field named ''%s''', fieldname);
        end
    end
    
    if verbose, fprintf('done. Saving.\n'); end
    % Save results to the file.
    save(save_file, 'results');
else
    if verbose, fprintf('Loading precomputed results...\t'); end
    contents = load(save_file);
    if verbose, fprintf('done.\n'); end
    results = contents.results;
end

varargout = results;
end