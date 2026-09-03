function test_mean_available()
%TEST_MEAN_AVAILABLE  mean() averages the channels the file has.
%
%   Exercises apply_channel_derivations directly with synthetic
%   (signal_header, signal_cell) inputs — no EDF file involved — so it
%   runs anywhere MATLAB does:
%
%       cd tests; test_mean_available
%
%   Covers: adaptive availability (a reference listing every cohort
%   spelling of a montage), bit-parity with the strict mean when all
%   arguments are present, argument all-or-nothing availability,
%   reference skip-with-warning when nothing is available, and the
%   per-Channels-entry skip for a channel using a skipped reference.

sig = @(v) v(:).';

% ---- 1. Partial availability: file has 2 of the 4 spellings ----------
sh = mk_sh({'C3-A2', 'O1-A2'});
c3 = sig([10 20 30]);  o1 = sig([4 6 8]);
sc = {c3, o1};
refs = {'M1 = mean($C3-A2$, $[C3-A2 - B]$, $O1-A2$, $[O1-A2 - B]$)'};
chans = {'C3 = $C3-A2$ - M1'};
[sh_out, sc_out] = apply_channel_derivations(sh, sc, chans, refs, false);
assert(numel(sh_out) == 1 && strcmp(sh_out(1).signal_labels, 'C3'));
expected = c3 - (c3 .* (1/2) + o1 .* (1/2));
assert(isequal(sc_out{1}, expected), 'partial-availability mean wrong');

% ---- 2. All arguments present: exact strict-mean arithmetic ----------
sh4 = mk_sh({'C3-A2', '[C3-A2 - B]', 'O1-A2', '[O1-A2 - B]'});
a = sig([1 2]); b = sig([3 4]); c = sig([5 6]); d = sig([7 8]);
sc4 = {a, b, c, d};
[sh_out, sc_out] = apply_channel_derivations(sh4, sc4, {'M = mean($C3-A2$, $[C3-A2 - B]$, $O1-A2$, $[O1-A2 - B]$)'}, {}, false);
assert(numel(sh_out) == 1);
strict = zeros(size(a));
for x = {a, b, c, d}
    strict = strict + (1/4) * x{1};   % the engine's own accumulation order
end
assert(isequal(sc_out{1}, strict), 'all-present mean must equal the strict 1/N form');

% ---- 3. An argument is all-or-nothing ---------------------------------
% mean(C3 - A2, O1): A2 missing, so the C3 - A2 argument must NOT
% contribute a half-referenced C3; the mean is just O1.
sh2 = mk_sh({'C3', 'O1'});
sc2 = {sig([10 10]), sig([4 4])};
[~, sc_out] = apply_channel_derivations(sh2, sc2, {'X = mean($C3$ - $A2$, $O1$)'}, {}, false);
assert(isequal(sc_out{1}, sig([4 4])), 'partially-present argument leaked into the mean');

% ---- 4. Nothing available: reference skipped, dependent channel skipped
w = warning('off', 'read_EDF:UnknownChannel');
cleanup = onCleanup(@() warning(w));
[sh_out, sc_out] = apply_channel_derivations(sh, sc, ...
    {'C4 = $C4-A1$ - M2', 'C3 = $C3-A2$ - M1'}, ...
    {'M2 = mean($C4-A1$, $[C4-A1 - B]$)', ...
     'M1 = mean($C3-A2$, $[C3-A2 - B]$)'}, false);
assert(numel(sh_out) == 1 && strcmp(sh_out(1).signal_labels, 'C3'), ...
    'unavailable reference must fail only its own channels');
assert(isequal(sc_out{1}, c3 - c3), 'single-spelling mean should be that channel');

% ---- 5. dry_run resolution matches the live path ----------------------
[sh_out, ~] = apply_channel_derivations(sh, cellfun(@(x) [], sc, 'UniformOutput', false), ...
    {'C4 = $C4-A1$ - M2', 'C3 = $C3-A2$ - M1'}, ...
    {'M2 = mean($C4-A1$, $[C4-A1 - B]$)', ...
     'M1 = mean($C3-A2$, $[C3-A2 - B]$)'}, false, true);
assert(numel(sh_out) == 1 && strcmp(sh_out(1).signal_labels, 'C3'), ...
    'dry_run must resolve the same channel set as a live run');

fprintf('test_mean_available: all assertions passed\n');
end


function sh = mk_sh(labels)
proto = struct('signal_labels', '', 'transducer_type', '', ...
    'physical_dimension', 'uV', 'physical_min', -100, 'physical_max', 100, ...
    'digital_min', -32768, 'digital_max', 32767, 'prefiltering', '', ...
    'samples_in_record', 100, 'sampling_frequency', 100);
sh = repmat(proto, 1, numel(labels));
for k = 1:numel(labels)
    sh(k).signal_labels = labels{k};
end
end
