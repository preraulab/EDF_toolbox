function profile_convert(in_dir, target_rate)
%PROFILE_CONVERT  Detailed phase-by-phase timing of the EDF convert pipeline.
%
%   profile_convert(in_dir, target_rate)
%
%   Walks through the same steps as convert_EDF but times each phase
%   independently:
%       read     : read_EDF
%       design   : resample filter coefficient design (first time per ratio)
%       run      : resample filter application
%       trim     : trim to record boundary
%       write    : write_EDF (quantize + compress + write)

if nargin < 2, target_rate = 100; end

toolbox = fileparts(mfilename('fullpath'));
addpath(genpath(toolbox));

files = dir(fullfile(in_dir, '*.edf'));
files = arrayfun(@(x) fullfile(x.folder, x.name), files, 'uni', false);
fprintf('=== profile_convert on %d files in %s ===\n', numel(files), in_dir);
fprintf('target_rate = %g Hz\n\n', target_rate);

filt_cache = containers.Map('KeyType','char','ValueType','any');
per_file = struct('file',{}, 'n_ch',{}, 'read',{}, 'design',{}, 'run',{}, ...
                  'trim',{}, 'write',{}, 'total',{});

for fi = 1:numel(files)
    f = files{fi};
    [~, base, ~] = fileparts(f);
    fprintf('--- [%d/%d] %s ---\n', fi, numel(files), base);

    % --- read ---
    t0 = tic;
    [header, signal_header, signal_cell, annotations] = read_EDF(f);
    t_read = toc(t0);

    % Locate annotation channel
    annot_idx = 0;
    for i = 1:numel(signal_header)
        if strcmp(strtrim(signal_header(i).signal_labels), 'EDF Annotations')
            annot_idx = i; break;
        end
    end
    record_duration = header.data_record_duration;
    n_signals = numel(signal_header);

    n_ch_resampled = 0;
    t_design = 0;
    t_run = 0;
    t_trim = 0;
    new_signal_cell = signal_cell;
    new_signal_header = signal_header;
    new_spr = round(target_rate * record_duration);

    for i = 1:n_signals
        if i == annot_idx, continue; end
        orig_rate = signal_header(i).sampling_frequency;
        if isempty(orig_rate) || ~isfinite(orig_rate) || orig_rate <= 0
            orig_rate = signal_header(i).samples_in_record / record_duration;
        end
        if abs(orig_rate - target_rate) < 1e-9, continue; end

        [P, Q] = rat(target_rate / orig_rate, 1e-6);
        key = sprintf('%d_%d', P, Q);
        x = signal_cell{i};

        if isKey(filt_cache, key)
            t1 = tic;
            y = resample(x(:), P, Q, filt_cache(key));
            t_run = t_run + toc(t1);
        else
            % Separate design and run by extracting filter from a tiny dummy
            % call (~design cost), then running the real signal with the
            % cached filter (~run cost).
            t1 = tic;
            [~, b] = resample(zeros(20,1), P, Q);
            t_design = t_design + toc(t1);
            filt_cache(key) = b;
            t1 = tic;
            y = resample(x(:), P, Q, b);
            t_run = t_run + toc(t1);
        end

        t1 = tic;
        y = y(:)';
        n_records = floor(numel(y) / new_spr);
        if numel(y) > n_records * new_spr
            y = y(1:n_records*new_spr);
        end
        new_signal_cell{i} = y;
        new_signal_header(i).samples_in_record  = new_spr;
        new_signal_header(i).sampling_frequency = target_rate;
        t_trim = t_trim + toc(t1);
        n_ch_resampled = n_ch_resampled + 1;
    end

    % Update header to match new sample counts
    min_records = inf;
    for i = 1:n_signals
        if i == annot_idx, continue; end
        r = floor(numel(new_signal_cell{i}) / new_signal_header(i).samples_in_record);
        if r < min_records, min_records = r; end
    end
    header.num_data_records = min_records;

    % --- write (zstd by suffix) ---
    out = fullfile(in_dir, ['_prof_' base '_' num2str(target_rate) 'Hz.edf.zst']);
    t1 = tic;
    write_EDF(out, header, new_signal_header, new_signal_cell, annotations, ...
        'AutoScale','recompute');
    t_write = toc(t1);
    if isfile(out), delete(out); end

    total = t_read + t_design + t_run + t_trim + t_write;
    fprintf('  channels resampled : %d\n', n_ch_resampled);
    fprintf('  read   : %6.3f s  (%5.1f%%)\n', t_read,   100*t_read/total);
    fprintf('  design : %6.3f s  (%5.1f%%)\n', t_design, 100*t_design/total);
    fprintf('  run    : %6.3f s  (%5.1f%%)\n', t_run,    100*t_run/total);
    fprintf('  trim   : %6.3f s  (%5.1f%%)\n', t_trim,   100*t_trim/total);
    fprintf('  write  : %6.3f s  (%5.1f%%)\n', t_write,  100*t_write/total);
    fprintf('  TOTAL  : %6.3f s\n\n', total);

    per_file(end+1) = struct('file', base, 'n_ch', n_ch_resampled, ...
        'read', t_read, 'design', t_design, 'run', t_run, ...
        'trim', t_trim, 'write', t_write, 'total', total); %#ok<AGROW>
end

% --- summary ---
T = struct2table(per_file);
fprintf('=== summary ===\n');
disp(T);
total_each = sum(T.read) + sum(T.design) + sum(T.run) + sum(T.trim) + sum(T.write);
fprintf('Aggregate over %d files: %.2f s\n', numel(files), total_each);
fprintf('  read   : %6.2f s  (%5.1f%%)\n', sum(T.read),   100*sum(T.read)/total_each);
fprintf('  design : %6.2f s  (%5.1f%%)  <- saved on cache hits\n', sum(T.design), 100*sum(T.design)/total_each);
fprintf('  run    : %6.2f s  (%5.1f%%)\n', sum(T.run),    100*sum(T.run)/total_each);
fprintf('  trim   : %6.2f s  (%5.1f%%)\n', sum(T.trim),   100*sum(T.trim)/total_each);
fprintf('  write  : %6.2f s  (%5.1f%%)\n', sum(T.write),  100*sum(T.write)/total_each);
end
