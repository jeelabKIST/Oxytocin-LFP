clear; clc;
addpath figure2/
addpath functions/

%% Load data
load('psd_time.mat');
load('drug.mat');
load ('cmap.mat');

%% Parameters
ch_names = {'mPFC', 'BLA', 'AC'};
ch_indices = [1, 2, 3];
target_ch_idx = [2, 1, 3];
target_ch_names = {'BLA', 'mPFC', 'AC'};
ch_names_fixed = {'AC', 'mPFC', 'BLA'};

band_list = {
    '\theta_{low}',  [4 8];
    '\theta_{high}', [8 12];
    '\beta_{low}',   [18 24];
    '\beta_{high}',  [24 32];
    '\gamma_{low}',  [35 50];
    '\gamma_{high}', [70 90]
    };
n_band = size(band_list, 1);

ot_list = find(drug == 1);
sal_list = find(drug == 0);

total_layers = n_band * 3;
dt_sample = 0.5;
color_matrix = [0.85 0.2 0.2; 0.55 0.2 0.2; 0.25 0.2 0.2];

%% Figure canvas
fig = figure(2); clf;
set(fig, 'Color', 'w', 'Position', [291 365 700 491]);

ax = subplot(3,3,1);
ax.Position = [0.0943 0.7093 0.2000 0.2157];
title('Surgery and LFP recordings');

%% Figure 2B
load('eeg_example.mat');

sb_x_len = 0.25;
sb_y_len = 0.1;
sb_text_y = '100 \muV';

toi = [3 5];
toi_idx = EEG.srate * 60 * toi(1) : EEG.srate * 60 * toi(2);
t = EEG.times(toi_idx);
doi = double(EEG.data(3, toi_idx));

[b_low, a_low] = butter(2, [8 12] / (EEG.srate/2), 'bandpass');
doi_low = filtfilt(b_low, a_low, doi);

[b_high, a_high] = butter(2, [30 50] / (EEG.srate/2), 'bandpass');
doi_high = filtfilt(b_high, a_high, doi);

subplot(3,3,2);
hold on;
plot(t, doi, 'k');
plot(t, doi_low - 0.16, 'k', 'LineWidth', 1);
plot(t, doi_high * 1.2 - 0.3, 'k', 'LineWidth', 1);
xlim([248.2 249.2]);
ylim([-0.4 0.2]);
axis off;

subplot(3,3,3);
hold on;

toi = [25 30];
toi_idx = EEG.srate * 60 * toi(1) : EEG.srate * 60 * toi(2);
t = EEG.times(toi_idx);
doi = double(EEG.data(3, toi_idx));

[b_low, a_low] = butter(2, [8 12] / (EEG.srate/2), 'bandpass');
doi_low = filtfilt(b_low, a_low, doi);

[b_high, a_high] = butter(2, [30 50] / (EEG.srate/2), 'bandpass');
doi_high = filtfilt(b_high, a_high, doi);

plot(t, doi, 'k');
plot(t, doi_low - 0.15, 'k', 'LineWidth', 1);
plot(t, doi_high - 0.3, 'k', 'LineWidth', 1);
xlim([1581.5 1582.5]);
ylim([-0.4 0.2]);
axis off;

xL = xlim;
text(xL(1)-0.045,  0.00, '- raw -', ...
    'HorizontalAlignment', 'right', 'FontSize', 11);
text(xL(1)-0.048, -0.15, '- \theta_{high}-', ...
    'HorizontalAlignment', 'right', 'Interpreter', 'tex', 'FontSize', 11);
text(xL(1)-0.055, -0.30, '- \gamma_{low}-', ...
    'HorizontalAlignment', 'right', 'Interpreter', 'tex', 'FontSize', 11);

x_start = 1582.5 - sb_x_len;
y_start = -0.4 + 0.05;

plot([x_start, x_start + sb_x_len], [y_start, y_start], '-k', 'LineWidth', 2);
plot([x_start + sb_x_len, x_start + sb_x_len], [y_start, y_start + sb_y_len], '-k', 'LineWidth', 2);

text(x_start + sb_x_len/2, y_start - 0.04, '250 ms', ...
    'HorizontalAlignment', 'center', 'FontSize', 12);
text(x_start + sb_x_len + 0.01, y_start + sb_y_len/2, sb_text_y, ...
    'HorizontalAlignment', 'left', 'FontSize', 12);

%% Figure 2C
subplot(3,3,4);
title('Power spectrogram');

% addpath('./main/stat');
load('psd_split.mat');

zpsd_split_set = (psd_split_set - m_psd_set) ./ s_psd_set;
spec_sal = zpsd_split_set(:,:,:,is_oxytocin == 0);
spec_oxy = zpsd_split_set(:,:,:,is_oxytocin == 1);

sz = size(spec_sal);
pvals = nan(sz(1:3));

pbar = progressbar(sz(1), "stat test");
for nr = 1:sz(1)
    for nc = 1:sz(2)
        for nch = 1:3
            x = squeeze(spec_sal(nr, nc, nch, :));
            y = squeeze(spec_oxy(nr, nc, nch, :));

            x = x(~isnan(x));
            y = y(~isnan(y));

            if isempty(x) || isempty(y)
                continue;
            end

            pvals(nr, nc, nch) = ranksum(x, y);
        end
    end
    pbar.update(nr);
end

alpha = 0.05;
for ch = 3
    im = median(spec_oxy(:,:,ch,:), 4, 'omitnan') - median(spec_sal(:,:,ch,:), 4, 'omitnan');

    alpha_mask = double(pvals(:,:,ch) < alpha);
    alpha_mask(alpha_mask == 0) = 0;
    alpha_mask(alpha_mask == 1) = 1;

    im_smooth = imgaussfilt(im, 1);
    alpha_smooth = imgaussfilt(alpha_mask, 1);

    interp_factor = 3;
    [tpsd_interp, fpsd_interp] = meshgrid( ...
        linspace(min(tpsd/60), max(tpsd/60), size(im, 2) * interp_factor), ...
        linspace(min(fpsd), max(fpsd), size(im, 1) * interp_factor));

    im_interp = interp2(tpsd/60, fpsd, im_smooth, tpsd_interp, fpsd_interp, 'linear');
    alpha_interp = interp2(tpsd/60, fpsd, alpha_smooth, tpsd_interp, fpsd_interp, 'linear');

    ax4 = subplot(3,3,4);
    ax4.Position = [0.0957 0.4299 0.1957 0.1933];

    imagesc(tpsd_interp(1,:), fpsd_interp(:,1), im_interp, 'AlphaData', alpha_interp);

    set(gca, 'LineWidth', 0.8, ...
        'XTick', [-19.5 0 20 40 59.5], ...
        'XTickLabel', {'-20','0','20','40','60'}, ...
        'YTick', [35 70]);

    colormap jet;
    cb = colorbar;
    ylabel(cb, 'z-scored power');
    cb.LineWidth = 0.8;
    cb.FontSize = 10;
    cb.Position = [0.3033 0.4297 0.0133 0.1935];
    cb.Ticks = [-0.6 0 0.6];

    clim([-0.6 0.6]);
    ylim([1 70]);
    title(sprintf('%s', ch_names{ch}));
    xlabel('Time (min)');
    ylabel('Freq (Hz)');
    xline(0, ':k', 'LineWidth', 1);
    box on;
    axis xy;
end

%% Figure 2E
subplot(3,3,[5 8]);
hold on;

dt_bin_z = 60;
pts_per_bin_z = dt_bin_z / dt_sample;

t_pre_vec = t_space_pre(:);
t_post_vec = t_space_post(:);

n_bin_pre = floor(length(t_pre_vec) / pts_per_bin_z);
n_bin_post = floor(length(t_post_vec) / pts_per_bin_z);

t_plot_z = [ ...
    linspace(t_pre_vec(1), t_pre_vec(end), n_bin_pre), ...
    linspace(t_post_vec(1), t_post_vec(end), n_bin_post)] / 60;

z_trace_mean = nan(length(t_plot_z), 3, n_band);

for c = 1:3
    for b = 1:n_band
        fidx = (freq >= band_list{b,2}(1)) & (freq <= band_list{b,2}(2));

        p_pre_raw  = squeeze(nanmean(d_pre(:, fidx, c, ot_list), 2));
        p_post_raw = squeeze(nanmean(d_post(:, fidx, c, ot_list), 2));

        p_pre_bin = squeeze(mean(reshape( ...
            p_pre_raw(end - n_bin_pre * pts_per_bin_z + 1:end, :), ...
            pts_per_bin_z, n_bin_pre, []), 1));

        p_post_bin = squeeze(mean(reshape( ...
            p_post_raw(1:n_bin_post * pts_per_bin_z, :), ...
            pts_per_bin_z, n_bin_post, []), 1));

        if size(p_pre_bin, 1) ~= n_bin_pre
            p_pre_bin = p_pre_bin';
        end
        if size(p_post_bin, 1) ~= n_bin_post
            p_post_bin = p_post_bin';
        end

        p_concat = [p_pre_bin; p_post_bin];
        mu_base = mean(p_concat(1:n_bin_pre, :), 1);
        sd_base = std(p_concat(1:n_bin_pre, :), 0, 1);

        z_trace_mean(:, c, b) = median((p_concat - mu_base) ./ sd_base, 2, 'omitnan');
    end
end

y_spacing = 1.2;
overlap_scale = 1.5;
xtick_pos = -20:20:60;

flipped_cmap = flipud(cmap{1}) * 0.95;
cmap_z = flipped_cmap;
cmap_rec = flipped_cmap;

idx_colors = round(linspace(1, size(cmap_z, 1), total_layers));
colors_z = cmap_z(idx_colors, :);
colors_rec = cmap_rec(idx_colors, :);

y_limit_max = total_layers * y_spacing;
y_min_limit = -y_spacing - 1.7;
y_max_limit = y_limit_max + 0.6;
tick_h = (y_max_limit - y_min_limit) * 0.01;

ax1 = subplot(3,3,[5 8]);
hold on;
ax1.Position = [0.4900 0.1100 0.2 0.5154];

tick_vals = [];
tick_labels = {};

for b = 1:n_band
    for k = 1:3
        c = target_ch_idx(k);
        layer_level = (b - 1) * 3 + k;
        y_base = (layer_level - 1) * y_spacing;

        smooth_trace = smoothdata(z_trace_mean(:, c, b), 'gaussian', 5);

        plot([-20 60], [y_base y_base], ':', 'Color', [0 0 0], 'LineWidth', 0.5);
        plot(t_plot_z, smooth_trace + y_base, 'Color', colors_z(layer_level, :), 'LineWidth', 1);

        tick_vals(end+1) = y_base; %#ok<SAGROW>
        tick_labels{end+1} = sprintf('%s, %s', target_ch_names{k}, band_list{b,1}); %#ok<SAGROW>
    end
end

xlabel('Time after injection (min)');
xlim([-10 60]);
ylim([y_min_limit y_max_limit]);

set(ax1, 'Box', 'off', 'XColor', 'w', 'YColor', 'w', ...
    'YTick', tick_vals, 'YTickLabel', tick_labels, 'FontSize', 10);

for x = xtick_pos
    line([x x], [y_min_limit, y_min_limit + tick_h], 'Color', 'k', 'LineWidth', 1);
end

ax1.XAxis.Label.Color = 'k';
ax1.YAxis.Label.Color = 'k';
ax1.XAxis.TickLabelColor = 'k';
ax1.YAxis.TickLabelColor = 'k';
ax1.XAxis.TickValues = xtick_pos;

%% Figure 2F
subplot(3,3,[6 9]);

dt_bin_rec = 120;
pts_per_bin_rec = dt_bin_rec / dt_sample;

t_all = [t_space_pre(:); t_space_post(:)];
n_bin_rec = floor(numel(t_all) / pts_per_bin_rec);
t_plot_rec = mean(reshape(t_all(1:n_bin_rec * pts_per_bin_rec), pts_per_bin_rec, n_bin_rec), 1) / 60;

pre_idx = t_plot_rec < 0;
inj_idx = find(pre_idx, 1, 'last');
t_inj = t_plot_rec(inj_idx);

t_recover_eachSess = nan(numel(ot_list), 3, n_band);
threshold = -1.28;

for c = 1:3
    for b = 1:n_band
        fidx = (freq >= band_list{b,2}(1)) & (freq <= band_list{b,2}(2));

        p_all = [ ...
            squeeze(nanmedian(d_pre(:, fidx, c, ot_list), 2)); ...
            squeeze(nanmedian(d_post(:, fidx, c, ot_list), 2))];

        p_binned = squeeze(nanmedian(reshape( ...
            p_all(1:n_bin_rec * pts_per_bin_rec, :), ...
            pts_per_bin_rec, n_bin_rec, []), 1));

        mu = nanmean(p_binned(pre_idx,:), 1);
        sd = nanstd(p_binned(pre_idx,:), 0, 1);
        z_b = (p_binned - mu) ./ sd;

        for s = 1:numel(ot_list)
            post_z = z_b(inj_idx+1:end, s);
            t_rel = t_plot_rec(inj_idx+1:end) - t_inj;
            supp_idx = find(post_z < threshold);

            if isempty(supp_idx)
                continue;
            end

            diff_idx = [0; diff(supp_idx) > 1];
            event_id = cumsum(diff_idx) + 1;

            t_rec = nan;

            for e = max(event_id):-1:1
                curr = supp_idx(event_id == e);
                if length(curr) >= 3
                    if curr(end) == length(post_z)
                        t_rec = 60;
                    else
                        t_rec = t_rel(curr(end) + 1);
                    end
                    break;
                end
            end

            t_recover_eachSess(s, c, b) = t_rec;
        end
    end
end

ax2 = subplot(3,3,[6 9]);
hold on;
ax2.Position = [0.7400 0.1100 0.2173 0.5154];

x_grid = linspace(0, 60, 200);
eps_margin = 1e-10;

for b = n_band:-1:1
    for k = 3:-1:1
        c = target_ch_idx(k);
        layer_level = (b - 1) * 3 + k;
        y_base = (layer_level - 1) * y_spacing;

        raw_vals = t_recover_eachSess(:, c, b);
        raw_vals = raw_vals(~isnan(raw_vals));

        if length(raw_vals) >= 2
            raw_vals(raw_vals < 0) = 0;
            raw_vals(raw_vals > 60) = 60;
            support_range = [0 - eps_margin, 60 + eps_margin];

            [f, ~] = ksdensity(raw_vals, x_grid, ...
                'Bandwidth', 4, ...
                'Support', support_range, ...
                'BoundaryCorrection', 'reflection');

            f = (f / (max(f) + 0.001)) * y_spacing * overlap_scale;

            fill([x_grid, 60, 0], [f + y_base, y_base, y_base], ...
                colors_rec(layer_level, :), ...
                'FaceAlpha', 0.7, ...
                'EdgeColor', [0.3 0.3 0.3], ...
                'LineWidth', 0.8);
        end

        plot([-20 60], [y_base y_base], ':', 'Color', [0 0 0], 'LineWidth', 0.5);
    end
end

xlabel('Time after injection (min)');
xlim([-10 60]);
ylim([y_min_limit y_max_limit]);

set(ax2, 'Box', 'off', 'XColor', 'w', 'YColor', 'w', ...
    'YTick', tick_vals, 'YTickLabel', [], 'FontSize', 10);

for x = xtick_pos
    line(ax2, [x x], [y_min_limit, y_min_limit + tick_h], 'Color', 'k', 'LineWidth', 1);
end

ax2.XAxis.Label.Color = 'k';
ax2.YAxis.Label.Color = 'k';
ax2.XAxis.TickLabelColor = 'k';
ax2.YAxis.TickLabelColor = 'w';
ax2.XAxis.TickValues = xtick_pos;

%% Figure 2D
subplot(3,3,7);

dA_oxy = cell(3, n_band);
dA_sal = cell(3, n_band);

df = [diff(freq(:)); 0];
idx_post = hb_findIdx([0 20] * 60, t_space_post);
idx_base = hb_findIdx([-5 0] * 60, t_space_pre);

for b = 1:n_band
    fidx = (freq >= band_list{b,2}(1)) & (freq <= band_list{b,2}(2));
    for ch = 1:3
        dA_oxy{ch,b} = band_delta_amp( ...
            squeeze(nanmean(d_post(idx_post,:,ch,ot_list),1)), ...
            squeeze(nanmean(d_pre(idx_base,:,ch,ot_list),1)), ...
            fidx, df) / 1000;

        dA_sal{ch,b} = band_delta_amp( ...
            squeeze(nanmean(d_post(idx_post,:,ch,sal_list),1)), ...
            squeeze(nanmean(d_pre(idx_base,:,ch,sal_list),1)), ...
            fidx, df) / 1000;
    end
end

ax3 = subplot(3,3,7);
hold on;
ax3.Position = [0.0957 0.1100 0.2477 0.2136];
title('Band specific changes');

bw = 0.22;
x_base = 1:6;
offsets = [-0.25, 0, 0.25];
h_bars = [];

y_min_limit = -0.8;
y_max_limit = 0.2;
tick_h = (y_max_limit - y_min_limit) * 0.03;

plot_order = [3 1 2];

for i = 1:3
    ch_idx = plot_order(i);
    for b = 1:6
        m = mean(dA_oxy{ch_idx,b}, 'omitnan');
        s = std(dA_oxy{ch_idx,b}, 'omitnan') / sqrt(nnz(isfinite(dA_oxy{ch_idx,b})));
        x_pos = x_base(b) + offsets(i);

        p = patch( ...
            [x_pos - bw/2, x_pos + bw/2, x_pos + bw/2, x_pos - bw/2], ...
            [0 0 m m], ...
            color_matrix(i,:), ...
            'FaceAlpha', 0.4, ...
            'EdgeColor', 'none');

        if b == 1
            h_bars(end+1) = p; %#ok<SAGROW>
        end

        plot([x_pos x_pos], [m - s m + s], 'Color', color_matrix(i,:), 'LineWidth', 1);

        [~, p_val] = ttest2(dA_oxy{ch_idx,b}, dA_sal{ch_idx,b});
        if p_val < (0.05 / 6)
            text(x_pos, max(0, m + s) + 0.02, '*', ...
                'Color', color_matrix(i,:), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 15);
        end
    end
end

yline(0, 'k:', 'LineWidth', 1);
ylabel('\Delta Amplitude (mV)');
ylim([y_min_limit y_max_limit]);
xlim([0.5 6.5]);

set(ax3, 'Box', 'off', 'XColor', 'w', 'YColor', 'w', ...
    'XTick', x_base, ...
    'XTickLabel', band_list(:,1), ...
    'TickLabelInterpreter', 'tex', ...
    'YTick', [-0.8 -0.4 0], ...
    'FontSize', 10);

for x = x_base
    line([x x], [y_min_limit, y_min_limit + tick_h], 'Color', 'k', 'LineWidth', 1);
end

y_ticks = [-0.8 -0.4 0];
for y = y_ticks
    line([0.5, 0.6], [y y], 'Color', 'k', 'LineWidth', 1);
end

ax3.XAxis.Label.Color = 'k';
ax3.YAxis.Label.Color = 'k';
ax3.XAxis.TickLabelColor = 'k';
ax3.YAxis.TickLabelColor = 'k';

lgd = legend(ax3, h_bars, ch_names_fixed, ...
    'Location', 'southeast', ...
    'Box', 'off', ...
    'FontSize', 8);
lgd.ItemTokenSize = [7, 3];

%% Local function
function dA = band_delta_amp(psd_post, psd_base, fidx, df)
    df_col = reshape(df(fidx), [], 1);
    P_post = squeeze(nansum(psd_post(fidx,:) .* df_col, 1));
    P_base = squeeze(nansum(psd_base(fidx,:) .* df_col, 1));
    dA = (sqrt(max(P_post,0)) - sqrt(max(P_base,0))) * 1000;
end