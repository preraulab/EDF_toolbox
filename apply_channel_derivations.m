function [sh_out, sc_out] = apply_channel_derivations(sh_in, sc_in, channels, references, verbose)
%APPLY_CHANNEL_DERIVATIONS  Build references and emit derived channels.
%
%   [SH_OUT, SC_OUT] = APPLY_CHANNEL_DERIVATIONS(SH_IN, SC_IN, CHANNELS,
%                                                REFERENCES, VERBOSE)
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
%   Errors: read_EDF:UnknownChannel, RateMismatch, BadMean,
%   RefCollision, ParseError. The 'read_EDF:' prefix is preserved
%   regardless of caller so existing error handlers keep working.

if nargin < 5, verbose = false; end
if nargin < 4, references = {}; end

% Step 1: pre-bake every reference into the loaded signal set.
[sh_aug, sc_aug] = compute_references(references, sh_in, sc_in, verbose);

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

    sig = zeros(size(sc_aug{leaf_idx(1)}));
    for j = 1:numel(spec.terms)
        sig = sig + spec.terms(j) * sc_aug{leaf_idx(j)};
    end

    new_sh = sh_aug(leaf_idx(1));
    new_sh.signal_labels = spec.label;
    new_sh.physical_min  = min(sig);
    new_sh.physical_max  = max(sig);

    sh_out(end+1) = new_sh; %#ok<AGROW>
    sc_out{end+1}  = sig;     %#ok<AGROW>

    % Chaining: aliased outputs become reusable names downstream.
    if ~isempty(user_alias)
        sh_aug(end+1)  = new_sh; %#ok<AGROW>
        sc_aug{end+1}  = sig;     %#ok<AGROW>
        aug_labels{end+1} = user_alias; %#ok<AGROW>
    end
end
end


function [sh_aug, sc_aug] = compute_references(refs, sh, sc, verbose)
%COMPUTE_REFERENCES  Pre-bake named References into the loaded signal set.
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

    wrapped = preprocess_spec(body, aug_labels);
    if verbose
        fprintf('apply_channel_derivations: ''%s = %s''  ->  ''%s = %s''\n', name, body, name, wrapped);
    end

    spec = parse_expr_string(wrapped, aug_labels);

    leaf_idx = zeros(1, numel(spec.leaves));
    for j = 1:numel(spec.leaves)
        leaf_idx(j) = find_label_idx(aug_labels, spec.leaves{j});
        if leaf_idx(j) == 0
            error('read_EDF:UnknownChannel', ...
                'Unknown channel ''%s'' in reference ''%s''.', spec.leaves{j}, s);
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
            'Reference ''%s'' has constituents at different sampling rates: %s', name, strjoin(parts, ', '));
    end

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
    ref_names_lower{end+1} = name_lower; %#ok<AGROW>
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

terms  = [];
leaves = {};

cur = 1;
[cur, sign] = parse_optional_sign(body, cur);

while true
    [cur, coefs, lvs] = parse_term(body, cur, all_labels, labels_lower, len_order, original);
    for j = 1:numel(lvs)
        terms(end+1)  = sign * coefs(j); %#ok<AGROW>
        leaves{end+1} = lvs{j};          %#ok<AGROW>
    end

    cur = skip_ws(body, cur);
    if cur > length(body)
        break
    end
    if body(cur) == '+'
        sign = +1; cur = cur + 1;
    elseif body(cur) == '-'
        sign = -1; cur = cur + 1;
    else
        error('read_EDF:ParseError', ...
            'Unexpected character ''%s'' at position %d in ''%s''.', body(cur), cur, original);
    end
end

if isempty(alias)
    out_label = strtrim(original);
else
    out_label = alias;
end

spec = struct('label', out_label, 'terms', terms, 'leaves', {leaves});
end


function [cur, coefs, leaves] = parse_term(body, cur, all_labels, labels_lower, len_order, original)
cur = skip_ws(body, cur);
if cur > length(body)
    error('read_EDF:ParseError', 'Unexpected end of expression in ''%s''.', original);
end

if cur + 3 <= length(body) && strcmpi(body(cur:cur+3), 'mean')
    peek = skip_ws(body, cur + 4);
    if peek <= length(body) && body(peek) == '('
        [cur, leaves] = parse_mean(body, peek + 1, all_labels, labels_lower, len_order, original);
        N = numel(leaves);
        if N < 2
            error('read_EDF:BadMean', ...
                'mean(...) needs at least 2 arguments in ''%s''.', original);
        end
        coefs = ones(1, N) / N;
        return
    end
end

[cur, label] = match_longest_label(body, cur, all_labels, labels_lower, len_order);
if isempty(label)
    error('read_EDF:UnknownChannel', ...
        'No matching channel at position %d in ''%s''.', cur, original);
end
coefs  = 1.0;
leaves = {label};
end


function [cur, leaves] = parse_mean(body, cur, all_labels, labels_lower, len_order, original)
leaves = {};
while true
    [cur, label] = match_longest_label(body, cur, all_labels, labels_lower, len_order);
    if isempty(label)
        error('read_EDF:UnknownChannel', ...
            'No matching channel inside mean() at position %d in ''%s''.', cur, original);
    end
    leaves{end+1} = label; %#ok<AGROW>
    cur = skip_ws(body, cur);
    if cur > length(body)
        error('read_EDF:ParseError', 'Unterminated mean(...) in ''%s''.', original);
    end
    if body(cur) == ','
        cur = cur + 1;
    elseif body(cur) == ')'
        cur = cur + 1;
        return
    else
        error('read_EDF:ParseError', ...
            'Expected '','' or '')'' inside mean() at position %d in ''%s''.', cur, original);
    end
end
end


function [cur, label] = match_longest_label(body, cur, all_labels, labels_lower, len_order)
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
