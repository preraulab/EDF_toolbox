function [sh_out, sc_out] = apply_channel_derivations(sh_in, sc_in, channels, references, verbose, dry_run)
%APPLY_CHANNEL_DERIVATIONS  Build references and emit derived channels.
%
%   [SH_OUT, SC_OUT] = APPLY_CHANNEL_DERIVATIONS(SH_IN, SC_IN, CHANNELS,
%                                                REFERENCES, VERBOSE, DRY_RUN)
%
%   DRY_RUN (logical, default false). When true, the function returns the
%   list of channel specs that WOULD have resolved (sh_out's labels) without
%   touching any signal data — only labels and sampling rates are needed.
%   sc_in may be an empty cell array (or a cellstr of {}/[] placeholders
%   length-aligned with sh_in); sc_out is returned as a cell of empties.
%   This is the resolution-only path used by the GUI's pre-flight scan to
%   answer "for this file's label set, which output names will I produce?"
%   without paying the cost of an EDF read. The same parser/RefCollision
%   semantics as the live path apply, so dry-run results match exactly what
%   a real run would emit.
%
%   Shared derivation pipeline used by both READ_EDF and WRITE_EDF so the
%   same expression syntax means the same thing on the read and write
%   sides. Inputs are an already-loaded (signal_header, signal_cell)
%   pair plus two cell arrays of expression strings. Outputs are the
%   transformed (signal_header, signal_cell) ready to be returned to
%   the user (read side) or fed into the writer (write side).
%
%   Pipeline:
%     1. compute_references evaluates 'References' in declared order,
%        appending each as a synthetic (signal_header, signal_cell)
%        entry. References are intermediates available to all
%        subsequent specs by name; they're dropped from the final
%        output unless the user also lists them in 'Channels'.
%
%     2. For each 'Channels' entry, in request order:
%          a. preprocess_spec wraps every real label occurrence
%             (longest-first) into '$LABEL$' tokens.
%          b. parse_expr_string parses the wrapped form into a flat
%             (terms, leaves) linear combination.
%          c. Leaves are resolved to indices, sampling rates are
%             checked, the signal is evaluated as
%             sig = sum_j coef_j * signal[leaf_j], and a new
%             signal_header is cloned from the first leaf.
%          d. If aliased ('OUT = expr'), the alias appends to the
%             namespace so subsequent entries can reference it.
%
%   When CHANNELS is empty: returns the original (non-synthetic)
%   channels untouched (References are dropped).
%
%   See HAS_DERIVED_SYNTAX for the cheap gate that decides whether the
%   pipeline is needed at all (callers route around it for plain label
%   subsets / legacy 'A-B' strings handled by the MEX backend).
%
%   mean(...) averages the AVAILABLE channels. Cohorts often name the
%   same montage differently across files, so a reference may list
%   every spelling — M = mean($C3-A2$, $[C3-A2 - B]$) — and each file
%   takes the mean over the arguments whose channels it actually has
%   (renormalized 1/M over the M available; error only when none are).
%   An argument counts only when every channel inside it is present.
%   When all N arguments are present the arithmetic is the exact 1/N
%   strict mean, unchanged. Quote each spelling in '$...$' so absent
%   ones parse as unresolved leaves rather than tripping ParseError.
%
%   Diagnostics. RateMismatch and BadMean always error.
%   UnknownChannel, ParseError, and RefCollision raised inside a
%   'Channels' entry are demoted to warnings — the offending entry is
%   skipped and the remaining channels still load. Inside a
%   'References' entry, structural problems (ParseError, RefCollision,
%   RateMismatch, BadMean) remain hard errors — they are wrong in every
%   file — but AVAILABILITY (UnknownChannel: the file simply lacks the
%   channels) skips that reference with a warning, and only the
%   'Channels' entries that use it then fail, per the usual per-entry
%   policy. The 'read_EDF:' prefix is preserved regardless of caller so
%   existing handlers keep working.

if nargin < 6, dry_run = false; end
if nargin < 5, verbose = false; end
if nargin < 4, references = {}; end

% Step 1: pre-bake every reference into the loaded signal set.
[sh_aug, sc_aug] = compute_references(references, sh_in, sc_in, verbose, dry_run);

n_real = numel(sh_in);

% Special case: References-only — emit the original channels untouched.
if isempty(channels)
    sh_out = sh_aug(1:n_real);
    sc_out = sc_aug(1:n_real);
    return
end

sh_out = sh_aug([]);
sc_out = {};

aug_labels = {sh_aug.signal_labels};

for k = 1:numel(channels)
    spec_str = strtrim(channels{k});

    try
        [user_alias, body] = split_alias(spec_str);

        % Preprocess only the BODY — the alias name is a fresh user
        % identifier and must not be searched for label substrings.
        wrapped_body = preprocess_spec(body, aug_labels);
        if isempty(user_alias)
            wrapped = wrapped_body;
        else
            wrapped = [user_alias ' = ' wrapped_body];
        end
        if verbose
            fprintf('apply_channel_derivations: ''%s''  ->  ''%s''\n', spec_str, wrapped);
        end

        spec = parse_expr_string(wrapped, aug_labels);

        % Restore the user's raw spec text as the output label when no alias.
        if isempty(user_alias)
            spec.label = spec_str;
        end

        % mean() over the channels this file has: drop the mean arguments
        % whose channels are absent, renormalize the rest. No-op when
        % everything (or nothing mean-shaped) is present.
        spec = resolve_mean_availability(spec, aug_labels, spec_str);

        leaf_idx = zeros(1, numel(spec.leaves));
        for j = 1:numel(spec.leaves)
            leaf_idx(j) = find_label_idx(aug_labels, spec.leaves{j});
            if leaf_idx(j) == 0
                error('read_EDF:UnknownChannel', ...
                    'Unknown channel ''%s'' in spec ''%s''.', spec.leaves{j}, spec_str);
            end
        end

        rates = zeros(1, numel(leaf_idx));
        for j = 1:numel(leaf_idx)
            rates(j) = sh_aug(leaf_idx(j)).sampling_frequency;
        end
        if any(rates ~= rates(1))
            parts = cell(1, numel(leaf_idx));
            for j = 1:numel(leaf_idx)
                parts{j} = sprintf('%s@%gHz', sh_aug(leaf_idx(j)).signal_labels, rates(j));
            end
            error('read_EDF:RateMismatch', ...
                'Sampling rates differ in spec ''%s'': %s', spec_str, strjoin(parts, ', '));
        end

        if ~isempty(user_alias) && find_label_idx(aug_labels, user_alias) > 0
            error('read_EDF:RefCollision', ...
                ['Channels alias ''%s'' collides with an existing channel ' ...
                 'or earlier output / reference.'], user_alias);
        end

        if dry_run
            % Resolution-only path: no signal arithmetic. Clone the
            % header, set the output label, and emit. Skip
            % physical_min/max (no signal to derive them from).
            new_sh = sh_aug(leaf_idx(1));
            new_sh.signal_labels = spec.label;
            sig = [];   %#ok<NASGU>  placeholder, not used
        else
            sig = zeros(size(sc_aug{leaf_idx(1)}));
            for j = 1:numel(spec.terms)
                sig = sig + spec.terms(j) * sc_aug{leaf_idx(j)};
            end

            new_sh = sh_aug(leaf_idx(1));
            new_sh.signal_labels = spec.label;
            new_sh.physical_min  = min(sig);
            new_sh.physical_max  = max(sig);
        end

        sh_out(end+1) = new_sh; %#ok<AGROW>
        if dry_run
            sc_out{end+1} = []; %#ok<AGROW>
        else
            sc_out{end+1}  = sig; %#ok<AGROW>
        end

        % Chaining: aliased outputs become reusable names downstream.
        % In dry_run mode we still need to extend aug_labels so later
        % specs that reference this alias resolve, but sc_aug just gets
        % a placeholder.
        if ~isempty(user_alias)
            sh_aug(end+1)  = new_sh; %#ok<AGROW>
            if dry_run
                sc_aug{end+1} = []; %#ok<AGROW>
            else
                sc_aug{end+1} = sig; %#ok<AGROW>
            end
            aug_labels{end+1} = user_alias; %#ok<AGROW>
        end
    catch ME
        % A per-entry resolution failure shouldn't take down the whole
        % load — warn and move on so the caller still gets every
        % channel that did resolve. Three cases are demoted to warnings:
        %   UnknownChannel : a leaf doesn't match any label.
        %   ParseError     : the parser tripped, typically because no
        %                    label could be wrapped — e.g. plain 'C3-A2'
        %                    against an EDF that only has '[C3-A2 - B]'.
        %   RefCollision   : an alias 'OUT = expr' tries to redefine a
        %                    name that already exists. Demoted so that
        %                    fallback patterns like {'A', 'A=B'} or
        %                    {'A=X', 'A=Y'} work the way you'd expect:
        %                    first match wins, later same-name entries
        %                    are quietly skipped.
        % Other errors (RateMismatch, BadMean) still propagate.
        if strcmp(ME.identifier, 'read_EDF:UnknownChannel') || ...
                strcmp(ME.identifier, 'read_EDF:ParseError') || ...
                strcmp(ME.identifier, 'read_EDF:RefCollision')
            % In dry_run mode the caller is asking "what would resolve?"
            % across many files — emitting a warning per dropped spec
            % per file would flood the command window. The drop itself
            % is reflected in sh_out (the spec doesn't appear there),
            % which is all the caller needs.
            if ~dry_run
                warning(ME.identifier, ...
                    'Skipping channel ''%s'': %s', spec_str, ME.message);
            end
            continue
        end
        rethrow(ME);
    end
end
end


function [sh_aug, sc_aug] = compute_references(refs, sh, sc, verbose, dry_run)
%COMPUTE_REFERENCES  Pre-bake named References into the loaded signal set.
if nargin < 5, dry_run = false; end
if nargin < 4, verbose = false; end

sh_aug = sh;
sc_aug = sc;

if isempty(refs)
    return
end

real_labels_lower = lower({sh.signal_labels});
ref_names_lower = {};

for k = 1:numel(refs)
    s = refs{k};

    [name, body] = split_alias(s);
    if isempty(name)
        error('read_EDF:ParseError', ...
            'Reference %d (''%s'') has no ''=''. References must be of the form ''NAME = expr''.', k, s);
    end
    if isempty(strtrim(body))
        error('read_EDF:ParseError', ...
            'Reference %d (''%s'') has empty right-hand side.', k, s);
    end

    name_lower = lower(name);
    if any(strcmp(real_labels_lower, name_lower)) || any(strcmp(ref_names_lower, name_lower))
        error('read_EDF:RefCollision', ...
            'Reference name ''%s'' collides with an existing channel or earlier reference.', name);
    end

    aug_labels = {sh_aug.signal_labels};

    % Availability is per file: a reference whose channels this file
    % simply doesn't have is SKIPPED (warning, not error), so cohorts
    % mixing montage naming schemes can declare one reference per
    % scheme. Only the output channels that use a skipped reference
    % fail — and those are already demoted per-entry by the caller.
    % Structural problems (ParseError, RefCollision, RateMismatch,
    % BadMean) stay hard errors: they are wrong in every file.
    try
        wrapped = preprocess_spec(body, aug_labels);
        if verbose
            fprintf('apply_channel_derivations: ''%s = %s''  ->  ''%s = %s''\n', name, body, name, wrapped);
        end

        spec = parse_expr_string(wrapped, aug_labels);

        % mean() over the channels this file has (see the helper).
        spec = resolve_mean_availability(spec, aug_labels, s);

        leaf_idx = zeros(1, numel(spec.leaves));
        for j = 1:numel(spec.leaves)
            leaf_idx(j) = find_label_idx(aug_labels, spec.leaves{j});
            if leaf_idx(j) == 0
                error('read_EDF:UnknownChannel', ...
                    'Unknown channel ''%s'' in reference ''%s''.', spec.leaves{j}, s);
            end
        end
    catch ME
        if strcmp(ME.identifier, 'read_EDF:UnknownChannel')
            if ~dry_run
                warning('read_EDF:UnknownChannel', ...
                    'Skipping reference ''%s'' (not available in this file): %s', name, ME.message);
            end
            continue
        end
        rethrow(ME);
    end

    rates = zeros(1, numel(leaf_idx));
    for j = 1:numel(leaf_idx)
        rates(j) = sh_aug(leaf_idx(j)).sampling_frequency;
    end
    if any(rates ~= rates(1))
        parts = cell(1, numel(leaf_idx));
        for j = 1:numel(leaf_idx)
            parts{j} = sprintf('%s@%gHz', sh_aug(leaf_idx(j)).signal_labels, rates(j));
        end
        error('read_EDF:RateMismatch', ...
            'Reference ''%s'' has constituents at different sampling rates: %s', name, strjoin(parts, ', '));
    end

    if dry_run
        new_sh = sh_aug(leaf_idx(1));
        new_sh.signal_labels = name;
        sh_aug(end+1) = new_sh; %#ok<AGROW>
        sc_aug{end+1} = []; %#ok<AGROW>
    else
        sig = zeros(size(sc_aug{leaf_idx(1)}));
        for j = 1:numel(spec.terms)
            sig = sig + spec.terms(j) * sc_aug{leaf_idx(j)};
        end

        new_sh = sh_aug(leaf_idx(1));
        new_sh.signal_labels = name;
        new_sh.physical_min  = min(sig);
        new_sh.physical_max  = max(sig);

        sh_aug(end+1) = new_sh; %#ok<AGROW>
        sc_aug{end+1}  = sig;     %#ok<AGROW>
    end
    ref_names_lower{end+1} = name_lower; %#ok<AGROW>
end
end


function spec = resolve_mean_availability(spec, aug_labels, where_str)
%RESOLVE_MEAN_AVAILABILITY  mean() averages the channels this file has.
%
%   Each mean(...) was parsed with its leaves stamped (grp/arg/grpn). An
%   argument is available only when EVERY leaf inside it resolves
%   against AUG_LABELS. Unavailable arguments are dropped and the kept
%   ones renormalized from the parse-time 1/N to 1/M (M = number
%   available), so a reference can list every cohort spelling of a
%   montage — mean($C3-A2$, $[C3-A2 - B]$) — and each file uses the
%   spellings it carries. When all N are available the spec is returned
%   untouched (bit-identical arithmetic to the historical strict mean).
%   A mean with NO available argument raises read_EDF:UnknownChannel,
%   which the callers demote per their usual per-entry policy. Leaves
%   outside any mean are left for the callers' own resolution errors.
if isempty(spec.grp) || all(spec.grp == 0)
    return
end

present = false(1, numel(spec.leaves));
for j = 1:numel(spec.leaves)
    present(j) = find_label_idx(aug_labels, spec.leaves{j}) > 0;
end

drop = false(1, numel(spec.leaves));
groups = unique(spec.grp(spec.grp > 0));
for gi = 1:numel(groups)
    g = groups(gi);
    in_g = spec.grp == g;
    args = unique(spec.arg(in_g));
    n_avail = 0;
    for ai = 1:numel(args)
        in_arg = in_g & spec.arg == args(ai);
        if all(present(in_arg))
            n_avail = n_avail + 1;
        else
            drop(in_arg) = true;
        end
    end
    if n_avail == 0
        error('read_EDF:UnknownChannel', ...
            'mean() in ''%s'': none of its channels are available (wanted any of ''%s'').', ...
            where_str, strjoin(spec.leaves(in_g), ''', '''));
    end
    N = spec.grpn(find(in_g, 1));
    if n_avail < N
        rescale = in_g & ~drop;
        spec.terms(rescale) = spec.terms(rescale) * (N / n_avail);
    end
end

if any(drop)
    spec.terms(drop)  = [];
    spec.leaves(drop) = [];
    spec.grp(drop)    = [];
    spec.arg(drop)    = [];
    spec.grpn(drop)   = [];
end
end


function wrapped = preprocess_spec(spec, all_labels)
%PREPROCESS_SPEC  Greedy longest-match wrap of real labels into $...$ tokens.
if isempty(spec) || isempty(all_labels)
    wrapped = spec;
    return
end

lengths = cellfun(@length, all_labels);
[~, order] = sort(lengths, 'descend');

spec_lower = lower(spec);
n = length(spec);
parts = {};
cur = 1;

while cur <= n
    ch = spec(cur);

    if ch == '$'
        rest = spec(cur+1:end);
        close_offsets = strfind(rest, '$');
        if isempty(close_offsets)
            parts{end+1} = spec(cur:end); %#ok<AGROW>
            cur = n + 1;
        else
            close_pos = cur + close_offsets(1);
            parts{end+1} = spec(cur:close_pos); %#ok<AGROW>
            cur = close_pos + 1;
        end
        continue
    end

    matched_label = '';
    matched_len = 0;
    for k = 1:numel(order)
        idx = order(k);
        L = all_labels{idx};
        Llen = length(L);
        if Llen == 0, continue; end
        if cur + Llen - 1 > n, continue; end
        if ~strcmp(spec_lower(cur:cur+Llen-1), lower(L)), continue; end
        % Word-boundary check on alphanumeric edges of both label and
        % surrounding spec, so a label like 'A' does NOT match inside
        % the word 'mean(', and a label 'AB' does not match inside
        % 'ABCD'. Pure-symbol labels (e.g. '[C3-A2]') are unaffected
        % because their first/last char is non-alphanumeric.
        L_first = L(1);
        L_last  = L(end);
        if is_alnum(L_first) && cur > 1 && is_alnum(spec(cur-1))
            continue
        end
        if is_alnum(L_last) && cur+Llen-1 < n && is_alnum(spec(cur+Llen))
            continue
        end
        matched_label = L;
        matched_len = Llen;
        break
    end

    if matched_len > 0
        parts{end+1} = ['$', matched_label, '$']; %#ok<AGROW>
        cur = cur + matched_len;
    else
        parts{end+1} = ch; %#ok<AGROW>
        cur = cur + 1;
    end
end

wrapped = strjoin(parts, '');
end


function [alias, body] = split_alias(s)
eq_pos = strfind(s, '=');
if isempty(eq_pos)
    alias = '';
    body  = s;
else
    alias = strtrim(s(1:eq_pos(1)-1));
    body  = strtrim(s(eq_pos(1)+1:end));
end
end


function idx = find_label_idx(all_labels, name)
target = lower(strtrim(name));
for i = 1:numel(all_labels)
    if strcmp(lower(strtrim(all_labels{i})), target)
        idx = i;
        return
    end
end
idx = 0;
end


function spec = parse_expr_string(s, all_labels)
% Parse a spec string into a flat (terms, leaves) linear combination.
%
% Grammar:
%   spec   := [name '='] expr
%   expr   := [sign] term (sign term)*
%   term   := factor (('*'|'/') factor)*
%   factor := ['+'|'-'] (number | '$' label '$' | mean '(' expr (',' expr)+ ')'
%                        | '(' expr ')')
%   sign   := '+' | '-'
%
% Coefficients fold at parse time: '(1/3)*C1 - 4*(C2-C3)/7' becomes
%   leaves = {C1, C2, C3}, terms = [1/3, -4/7, 4/7]. The engine then
% computes  sig = sum_j terms(j) * signal[leaves(j)].
%
% The grammar is restricted to LINEAR combinations of signals: a term
% may contain at most one signal-valued factor. Scalar*scalar and
% scalar*(signal-expr) and (signal-expr)/scalar fold as expected;
% scalar/(signal-expr), scalar+(signal-expr), or two signal-valued
% factors multiplied/divided are rejected as non-linear.

original = s;
[alias, body] = split_alias(s);
body = strtrim(body);

if isempty(body)
    error('read_EDF:ParseError', 'Empty expression: ''%s''', original);
end

labels_lower = cell(1, numel(all_labels));
for k = 1:numel(all_labels)
    labels_lower{k} = lower(strtrim(all_labels{k}));
end
[~, len_order] = sort(cellfun(@length, labels_lower), 'descend');

cur = 1;
[result, cur] = parse_expr(body, cur, all_labels, labels_lower, len_order, original);

cur = skip_ws(body, cur);
if cur <= length(body)
    error('read_EDF:ParseError', ...
        'Unexpected character ''%s'' at position %d in ''%s''.', ...
        body(cur), cur, original);
end

if result.is_scalar
    error('read_EDF:ParseError', ...
        'Expression ''%s'' has no signal — pure scalar.', original);
end

if isempty(alias)
    out_label = strtrim(original);
else
    out_label = alias;
end

spec = struct('label', out_label, ...
    'terms', result.coefs, 'leaves', {result.leaves}, ...
    'grp', result.grp, 'arg', result.arg, 'grpn', result.grpn);
end


function [r, cur] = parse_expr(body, cur, all_labels, labels_lower, len_order, original)
[cur, sgn] = parse_optional_sign(body, cur);
[r, cur]   = parse_term(body, cur, all_labels, labels_lower, len_order, original);
r = scale_result(r, sgn);

while true
    cur = skip_ws(body, cur);
    if cur > length(body) || (body(cur) ~= '+' && body(cur) ~= '-')
        return
    end
    if body(cur) == '+', sgn = 1; else, sgn = -1; end
    cur = cur + 1;
    [t, cur] = parse_term(body, cur, all_labels, labels_lower, len_order, original);
    t = scale_result(t, sgn);
    r = combine_add(r, t, original);
end
end


function [r, cur] = parse_term(body, cur, all_labels, labels_lower, len_order, original)
[r, cur] = parse_factor(body, cur, all_labels, labels_lower, len_order, original);
while true
    cur = skip_ws(body, cur);
    if cur > length(body) || (body(cur) ~= '*' && body(cur) ~= '/')
        return
    end
    op  = body(cur);
    cur = cur + 1;
    [f, cur] = parse_factor(body, cur, all_labels, labels_lower, len_order, original);
    r = combine_muldiv(r, op, f, original);
end
end


function [r, cur] = parse_factor(body, cur, all_labels, labels_lower, len_order, original)
cur = skip_ws(body, cur);
[cur, sgn] = parse_optional_sign(body, cur);

cur = skip_ws(body, cur);
if cur > length(body)
    error('read_EDF:ParseError', 'Unexpected end of expression in ''%s''.', original);
end

ch = body(cur);

if ch == '('
    [inner, cur] = parse_expr(body, cur + 1, all_labels, labels_lower, len_order, original);
    cur = skip_ws(body, cur);
    if cur > length(body) || body(cur) ~= ')'
        error('read_EDF:ParseError', 'Missing '')'' in ''%s''.', original);
    end
    cur = cur + 1;
    r = inner;

elseif ch == '$'
    [cur, label, quoted_raw] = match_longest_label(body, cur, all_labels, labels_lower, len_order);
    if isempty(label)
        if isempty(quoted_raw)
            error('read_EDF:UnknownChannel', ...
                'No matching channel at position %d in ''%s''.', cur, original);
        end
        % A well-formed '$...$' token naming a channel this file does not
        % have: keep it as an unresolved leaf instead of failing the
        % parse. If it sits inside a mean(...), the availability pass
        % drops that argument; anywhere else, leaf resolution raises the
        % same read_EDF:UnknownChannel it always did.
        r = vec_factor(1.0, quoted_raw);
    else
        r = vec_factor(1.0, label);
    end

elseif (ch >= '0' && ch <= '9') || ch == '.'
    [val, cur] = parse_number(body, cur, original);
    r = scalar_factor(val);

elseif cur + 3 <= length(body) && strcmpi(body(cur:cur+3), 'mean')
    peek = skip_ws(body, cur + 4);
    if peek <= length(body) && body(peek) == '('
        [r, cur] = parse_mean(body, peek + 1, all_labels, labels_lower, len_order, original);
    else
        % 'mean' not followed by '(' — try as a label match (longest-prefix
        % label that happens to start with 'mean'). The preprocessor
        % normally wraps real labels in $...$, so this path is unusual.
        [cur, label] = match_longest_label(body, cur, all_labels, labels_lower, len_order);
        if isempty(label)
            error('read_EDF:ParseError', ...
                'Expected ''('' after ''mean'' at position %d in ''%s''.', cur, original);
        end
        r = vec_factor(1.0, label);
    end

else
    % Last-ditch: try a bare longest-prefix label match (defensive — the
    % preprocessor should already have wrapped real labels in $...$).
    [cur2, label] = match_longest_label(body, cur, all_labels, labels_lower, len_order);
    if ~isempty(label)
        cur = cur2;
        r = vec_factor(1.0, label);
    else
        error('read_EDF:ParseError', ...
            'Unexpected character ''%s'' at position %d in ''%s''.', ch, cur, original);
    end
end

r = scale_result(r, sgn);
end


function [r, cur] = parse_mean(body, cur, all_labels, labels_lower, len_order, original)
% Caller has already consumed 'mean' and '('. Parse (expr (',' expr)+) ')'.
args = {};
while true
    [a, cur] = parse_expr(body, cur, all_labels, labels_lower, len_order, original);
    args{end+1} = a; %#ok<AGROW>
    cur = skip_ws(body, cur);
    if cur > length(body)
        error('read_EDF:ParseError', 'Unterminated mean(...) in ''%s''.', original);
    end
    if body(cur) == ','
        cur = cur + 1;
    elseif body(cur) == ')'
        cur = cur + 1;
        break
    else
        error('read_EDF:ParseError', ...
            'Expected '','' or '')'' inside mean() at position %d in ''%s''.', cur, original);
    end
end

N = numel(args);
if N < 2
    error('read_EDF:BadMean', ...
        'mean(...) needs at least 2 arguments in ''%s''.', original);
end

% mean = (1/N) * sum_i args_i, evaluated over the arguments the file
% actually has: every leaf is stamped with (group, argument, N) so
% resolve_mean_availability can drop the arguments whose channels are
% absent and renormalize the rest to 1/M. With all N available nothing
% is touched, so fully-covered files keep the exact 1/N arithmetic.
% The group id is this mean's character position in the spec — unique
% within one parse. Stamping OVERWRITES any marks from a nested inner
% mean: an inner mean is all-or-nothing within its enclosing argument
% (only the outermost mean adapts).
g = cur;
r = empty_vec();
for k = 1:N
    a = args{k};
    if a.is_scalar
        error('read_EDF:BadMean', ...
            'mean() argument %d is a pure scalar in ''%s''.', k, original);
    end
    a = scale_result(a, 1/N);
    a.grp(:)  = g;
    a.arg(:)  = k;
    a.grpn(:) = N;
    r = combine_add(r, a, original);
end
end


function [val, cur] = parse_number(body, cur, original)
% Match decimal number with optional fraction and exponent.
n = length(body);
start = cur;
while cur <= n && body(cur) >= '0' && body(cur) <= '9'
    cur = cur + 1;
end
if cur <= n && body(cur) == '.'
    cur = cur + 1;
    while cur <= n && body(cur) >= '0' && body(cur) <= '9'
        cur = cur + 1;
    end
end
if cur <= n && (body(cur) == 'e' || body(cur) == 'E')
    cur = cur + 1;
    if cur <= n && (body(cur) == '+' || body(cur) == '-')
        cur = cur + 1;
    end
    while cur <= n && body(cur) >= '0' && body(cur) <= '9'
        cur = cur + 1;
    end
end
str = body(start:cur-1);
val = str2double(str);
if isnan(val)
    error('read_EDF:ParseError', 'Bad number ''%s'' in ''%s''.', str, original);
end
end


% Result structs carry three arrays parallel to coefs/leaves recording
% mean() membership: grp (mean-group id, 0 = not inside a mean), arg
% (argument index within the group), grpn (the group's total argument
% count N). resolve_mean_availability uses them to drop the arguments a
% file doesn't have and renormalize the rest — see that function.
function r = scalar_factor(v)
r = struct('is_scalar', true, 'scalar', v, 'coefs', [], 'leaves', {{}}, ...
    'grp', [], 'arg', [], 'grpn', []);
end


function r = vec_factor(coef, label)
r = struct('is_scalar', false, 'scalar', 0, 'coefs', coef, 'leaves', {{label}}, ...
    'grp', 0, 'arg', 0, 'grpn', 0);
end


function r = empty_vec()
r = struct('is_scalar', false, 'scalar', 0, 'coefs', [], 'leaves', {{}}, ...
    'grp', [], 'arg', [], 'grpn', []);
end


function r = scale_result(r, s)
if s == 1, return; end
if r.is_scalar
    r.scalar = r.scalar * s;
else
    r.coefs = r.coefs * s;
end
end


function r = combine_add(a, b, original)
if a.is_scalar && b.is_scalar
    r = scalar_factor(a.scalar + b.scalar);
elseif a.is_scalar || b.is_scalar
    error('read_EDF:ParseError', ...
        'Cannot add a scalar to a signal expression in ''%s'' (DC offsets not supported).', original);
else
    r = struct('is_scalar', false, 'scalar', 0, ...
        'coefs', [a.coefs, b.coefs], ...
        'leaves', {[a.leaves, b.leaves]}, ...
        'grp',  [a.grp,  b.grp], ...
        'arg',  [a.arg,  b.arg], ...
        'grpn', [a.grpn, b.grpn]);
end
end


function r = combine_muldiv(a, op, b, original)
if a.is_scalar && b.is_scalar
    if op == '*'
        r = scalar_factor(a.scalar * b.scalar);
    else
        if b.scalar == 0
            error('read_EDF:ParseError', 'Division by zero in ''%s''.', original);
        end
        r = scalar_factor(a.scalar / b.scalar);
    end
elseif a.is_scalar && ~b.is_scalar
    if op == '*'
        r = scale_result(b, a.scalar);
    else
        error('read_EDF:ParseError', ...
            'Cannot divide a scalar by a signal expression in ''%s'' (not linear).', original);
    end
elseif ~a.is_scalar && b.is_scalar
    if op == '*'
        r = scale_result(a, b.scalar);
    else
        if b.scalar == 0
            error('read_EDF:ParseError', 'Division by zero in ''%s''.', original);
        end
        r = scale_result(a, 1.0 / b.scalar);
    end
else
    if op == '*'
        verb = 'multiply';
    else
        verb = 'divide';
    end
    error('read_EDF:ParseError', ...
        'Cannot %s two signal expressions in ''%s'' (not linear).', verb, original);
end
end


function [cur, label, quoted_raw] = match_longest_label(body, cur, all_labels, labels_lower, len_order)
% Third output: the raw text of a well-formed '$...$' token when it did
% not match any label ('' otherwise) — lets the caller keep an unknown
% quoted channel as an unresolved leaf rather than failing the parse.
quoted_raw = '';
cur = skip_ws(body, cur);
if cur > length(body)
    label = '';
    return
end

if body(cur) == '$'
    end_pos = strfind(body(cur+1:end), '$');
    if isempty(end_pos)
        error('read_EDF:ParseError', ...
            'Unterminated $-quoted label at position %d in ''%s''.', cur, body);
    end
    quoted = strtrim(body(cur+1 : cur+end_pos(1)-1));
    quoted_lower = lower(quoted);
    label = '';
    for k = 1:numel(all_labels)
        if strcmp(lower(strtrim(all_labels{k})), quoted_lower) %#ok<STCI>
            label = all_labels{k};
            break
        end
    end
    if isempty(label)
        quoted_raw = quoted;
    end
    cur = cur + end_pos(1) + 1;
    return
end

rest_lower = lower(body(cur:end));
for k = 1:numel(len_order)
    idx = len_order(k);
    L = labels_lower{idx};
    if isempty(L), continue; end
    if length(rest_lower) >= length(L) && strcmp(rest_lower(1:length(L)), L)
        label = all_labels{idx};
        cur   = cur + length(L);
        return
    end
end
label = '';
end


function cur = skip_ws(s, cur)
while cur <= length(s) && (s(cur) == ' ' || s(cur) == sprintf('\t'))
    cur = cur + 1;
end
end


function [cur, sign] = parse_optional_sign(s, cur)
cur = skip_ws(s, cur);
if cur <= length(s) && s(cur) == '-'
    sign = -1; cur = cur + 1;
elseif cur <= length(s) && s(cur) == '+'
    sign = +1; cur = cur + 1;
else
    sign = +1;
end
end


function tf = is_alnum(c)
tf = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
end
