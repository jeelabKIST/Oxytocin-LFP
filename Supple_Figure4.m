%% =========================================
% Old-style time-resolved 5-panel plot
% aligned to injection end
%
% Panels:
% 1) mPFC low-theta power
% 2) AC low-theta power
% 3) mPFC-AC low-theta PLV
% 4) mPFC -> AC low-theta GC
% 5) AC -> mPFC low-theta GC
%
% Features:
% - same valid-session logic as old GC code
% - manual outlier removal
% - successive 5-min windows
% - significance shown as one merged bar + one star
% - inset bar plot for 0-20 min average
%% =========================================

clearvars
close all; clc

%% -----------------------------------------
% 1) load data
% -----------------------------------------
load('data/processed/gc_cleaned.mat')      % gc, meta_info
load('data/processed/psd_time_250511d.mat') % psd
load('./data/processed/plv_0401b.mat')      % plv

inj_time   = meta_info.inj_time;
drug       = meta_info.drug;
skip_ratio = meta_info.skip_ratio;

%% -----------------------------------------
% 2) valid sessions: same spirit as old code
% -----------------------------------------
isEmptyGC  = cellfun(@isempty, gc(:,1));
isEmptyPSD = cellfun(@isempty, psd(:,1));
isEmptyPLV = cellfun(@isempty, plv(:,1));

valid_fileIdx = find((skip_ratio <= 10) & ~isEmptyGC & ~isEmptyPSD & ~isEmptyPLV);

drug_valid     = drug(valid_fileIdx);
inj_time_valid = inj_time(valid_fileIdx,:);

%% -----------------------------------------
% 3) settings
% -----------------------------------------
% channels: mPFC-BLA-AC
ch_mPFC = 1;
ch_AC   = 3;

% pair index: 1=mPFC-BLA, 2=mPFC-AC, 3=BLA-AC
pair_mPFC_AC = 2;

% GC direction assumption from your old code:
% dir 1 = mPFC -> AC
% dir 2 = AC -> mPFC
% If opposite, swap labels below.
dir_mPFC_to_AC = 1;
dir_AC_to_mPFC = 2;

% low-theta band
theta_range = [4 8];

% frequency vectors
freq_psd = 1:121;
freq_plv = 1:100;
freq_gc  = 1:512;

f_psd = freq_psd >= theta_range(1) & freq_psd <= theta_range(2);
f_plv = freq_plv >= theta_range(1) & freq_plv <= theta_range(2);
f_gc  = freq_gc  >= theta_range(1) & freq_gc  <= theta_range(2);

% old outlier removal
sal_outliers = [4, 13, 14, 16, 18, 19];
oxy_outliers = [ 2, 16, 16:18];

sal_idx = find(drug_valid == 0);
oxy_idx = find(drug_valid == 1);

sal_idx(sal_outliers) = [];
oxy_idx(oxy_outliers) = [];

% successive 5-min windows aligned to injection end
toi_edges = -20*60 : 5*60 : 60*60;   % sec relative to injection end
n_bin = numel(toi_edges) - 1;
time_axis = (toi_edges(1:end-1) + toi_edges(2:end)) / 2 / 60;

%% -----------------------------------------
% 4) storage
% -----------------------------------------
% session x bin
pow_mPFC_all = nan(length(valid_fileIdx), n_bin);
pow_AC_all   = nan(length(valid_fileIdx), n_bin);
plv_all      = nan(length(valid_fileIdx), n_bin);
gc_fw_all    = nan(length(valid_fileIdx), n_bin); % mPFC -> AC
gc_bw_all    = nan(length(valid_fileIdx), n_bin); % AC -> mPFC

%% -----------------------------------------
% 5) compute time-resolved values
% -----------------------------------------
for newIdx = 1:length(valid_fileIdx)
    fileIdx = valid_fileIdx(newIdx);

    gc_data  = real(gc{fileIdx,1});
    t_gc     = gc{fileIdx,2};

    psd_data = psd{fileIdx,1};
    t_psd    = psd{fileIdx,2};

    plv_data = plv{fileIdx,1};
    t_plv    = plv{fileIdx,2};

    if isempty(gc_data) || isempty(t_gc) || isempty(psd_data) || isempty(t_psd) || isempty(plv_data) || isempty(t_plv)
        continue
    end

    if inj_time_valid(newIdx,1) <= 0
        continue
    end

    % align to injection end using gc time base
    inj_start_end_gc = hb_findIdx(inj_time_valid(newIdx,:), t_gc);
    inj_end_time = t_gc(inj_start_end_gc(end));

    for b = 1:n_bin
        win_start_t = inj_end_time + toi_edges(b);
        win_end_t   = inj_end_time + toi_edges(b+1);

        % --- PSD window
        idx_psd = t_psd >= win_start_t & t_psd <= win_end_t;
        if any(idx_psd)
            tmp = psd_data(idx_psd, f_psd, ch_mPFC);
            tmp(~isfinite(tmp)) = NaN;
            pow_mPFC_all(newIdx,b) = mean(tmp(:), 'omitnan');

            tmp = psd_data(idx_psd, f_psd, ch_AC);
            tmp(~isfinite(tmp)) = NaN;
            pow_AC_all(newIdx,b) = mean(tmp(:), 'omitnan');
        end

        % --- PLV window
        idx_plv = t_plv >= win_start_t & t_plv <= win_end_t;
        if any(idx_plv)
            tmp = plv_data(idx_plv, f_plv, pair_mPFC_AC);
            tmp(~isfinite(tmp)) = NaN;
            plv_all(newIdx,b) = mean(tmp(:), 'omitnan');
        end

        % --- GC window
        idx_gc = t_gc >= win_start_t & t_gc <= win_end_t;
        if any(idx_gc)
            indiv_gc = squeeze(nanmedian(gc_data(idx_gc,:,:,:), 1)); % freq x pair x dir

            tmp = indiv_gc(:, pair_mPFC_AC, dir_mPFC_to_AC);
            tmp(~isfinite(tmp) | tmp == 0) = NaN;
            gc_fw_all(newIdx,b) = mean(tmp(f_gc), 'omitnan');

            tmp = indiv_gc(:, pair_mPFC_AC, dir_AC_to_mPFC);
            tmp(~isfinite(tmp) | tmp == 0) = NaN;
            gc_bw_all(newIdx,b) = mean(tmp(f_gc), 'omitnan');
        end
    end
end

%% -----------------------------------------
% 6) group matrices
% -----------------------------------------
data_all = {pow_mPFC_all, pow_AC_all, plv_all, gc_fw_all, gc_bw_all};

valid_sal_panel = cell(1,5);
valid_oxy_panel = cell(1,5);

for p = 1:5
    sal_mat = data_all{p}(sal_idx,:);
    oxy_mat = data_all{p}(oxy_idx,:);

    % validity filtering
    valid_sal = all(~isnan(sal_mat), 2) & any(sal_mat ~= 0, 2);
    valid_oxy = all(~isnan(oxy_mat), 2) & any(oxy_mat ~= 0, 2);

    valid_sal_panel{p} = sal_mat(valid_sal,:);
    valid_oxy_panel{p} = oxy_mat(valid_oxy,:);
end

%% -----------------------------------------
% 7) manual point-level cleanup
% -----------------------------------------
% 40 min 부근 bin
% bad_bin = find(time_axis >= 37.5 & time_axis <= 42.5);

% 사용 예시:
% 패널 번호:
% 1=mPFC power, 2=AC power, 3=PLV, 4=mPFC->AC GC, 5=AC->mPFC GC

% A. 특정 패널의 해당 bin 전체 제거
% for p = [4 5]
%     valid_oxy_panel{p}(:, bad_bin) = NaN;
%     valid_sal_panel{p}(:, bad_bin) = NaN;
% end

% B. 특정 세션만 제거
% suspicious_oxy_sessions = [3 7];
% suspicious_sal_sessions = [5];
% p = 5;
% valid_oxy_panel{p}(suspicious_oxy_sessions, bad_bin) = NaN;
% valid_sal_panel{p}(suspicious_sal_sessions, bad_bin) = NaN;

%% -----------------------------------------
% 8) plot
% -----------------------------------------
fig = figure(1);clf
set(fig, 'Position',[100 100 1400 300],'Color','w');

panel_names = { ...
    'mPFC low-\theta power', ...
    'AC low-\theta power', ...
    'mPFC-AC PLV (\theta_{low})', ...
    'mPFC \rightarrow AC GC (\theta_{low})', ...
    'AC \rightarrow mPFC GC (\theta_{low})'};

% ylabels = { ...
%     'Power (a.u.)', ...
%     'Power (a.u.)', ...
%     'PLV', ...
%     'Low-\theta GC (a.u.)', ...
%     'Low-\theta GC (a.u.)'};

ylabels = { ...
    '\Delta Power (a.u.)', ...
    '\Delta Power (a.u.)', ...
    '\Delta PLV', ...
    '\Delta Low-\theta GC (a.u.)', ...
    '\Delta GC (a.u.)'};

col_sal = [0.3 0.3 0.3];
col_oxy = [0.85 0.2 0.2];

for p = 1:5
    sub_ax = subplot(1,5,p); hold on; axis square

    sal_mat = valid_sal_panel{p};
    oxy_mat = valid_oxy_panel{p};


    % ---------- baseline normalization using pre-injection bins ----------
    base_idx = time_axis < 0;   % pre-injection bins

    sal_base = mean(sal_mat(:, base_idx), 2, 'omitnan');   % session x 1
    oxy_base = mean(oxy_mat(:, base_idx), 2, 'omitnan');

    % subtract session-specific baseline
    sal_mat = sal_mat - sal_base;
    oxy_mat = oxy_mat - oxy_base;
        % ---------- baseline normalization using pre-injection bins ----------


    mean_sal = mean(sal_mat, 1, 'omitnan');
    sem_sal  = std(sal_mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(sal_mat),1));

    mean_oxy = mean(oxy_mat, 1, 'omitnan');
    sem_oxy  = std(oxy_mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(oxy_mat),1));

    % shaded lines
    fill([time_axis fliplr(time_axis)], ...
         [mean_sal+sem_sal fliplr(mean_sal-sem_sal)], ...
         [0.6 0.6 0.6], 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(time_axis, mean_sal, '-', 'Color',col_sal, 'LineWidth',1.5);

    fill([time_axis fliplr(time_axis)], ...
         [mean_oxy+sem_oxy fliplr(mean_oxy-sem_oxy)], ...
         [1 0.3 0.3], 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(time_axis, mean_oxy, '-', 'Color',col_oxy, 'LineWidth',1.5);

    xline(0, 'k--');

    % ---------- significance as merged bar + one star ----------
    sig_bins = false(1,n_bin);

    for b = 1:n_bin
        x1 = oxy_mat(:,b);
        x2 = sal_mat(:,b);

        x1 = x1(isfinite(x1));
        x2 = x2(isfinite(x2));

        if numel(x1) < 3 || numel(x2) < 3
            continue
        end

        [~, pval] = ttest2(x1, x2);
        if pval < 0.05
            sig_bins(b) = true;
        end
    end

    d = diff([0 sig_bins 0]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;

    y_lim = ylim;
    y_sig = y_lim(2) - 0.06*(y_lim(2)-y_lim(1));

    for k = 1:length(start_idx)
        b1 = start_idx(k);
        b2 = end_idx(k);

        x_start = toi_edges(b1)   / 60;
        x_end   = toi_edges(b2+1) / 60;

        x1 = mean(oxy_mat(:,b1:b2), 2, 'omitnan');
        x2 = mean(sal_mat(:,b1:b2), 2, 'omitnan');

        sig_col = 'k';
        if mean(x1,'omitnan') > mean(x2,'omitnan')
            sig_col = 'r';
        end

        line([x_start x_end], [y_sig y_sig], 'Color', sig_col, 'LineWidth', 2);
        text((x_start+x_end)/2, y_sig, '*', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',11, ...
            'Color',sig_col);
    end

    xlabel('Time (min)');
    ylabel(ylabels{p});
    title(panel_names{p}, 'FontWeight','normal');

    set(gca, 'Box','off', 'TickDir','out', 'LineWidth',1, ...
        'FontSize',10, 'FontName','Arial');

    % ---------- inset bar: 5-20 min average ----------
    idx_avg = time_axis >= 0 & time_axis <= 20;

    val_sal = mean(sal_mat(:, idx_avg), 2, 'omitnan');
    val_oxy = mean(oxy_mat(:, idx_avg), 2, 'omitnan');

    val_sal = val_sal(isfinite(val_sal));
    val_oxy = val_oxy(isfinite(val_oxy));

    pos = get(sub_ax, 'Position');
    inset_pos = [pos(1)+pos(3)*0.56, pos(2)+pos(4)*0.55, pos(3)*0.32, pos(4)*0.28];
    axes('Position', inset_pos); hold on;

    m_bars = [mean(val_oxy,'omitnan'), mean(val_sal,'omitnan')];
    s_bars = [std(val_oxy,'omitnan')/sqrt(numel(val_oxy)), ...
              std(val_sal,'omitnan')/sqrt(numel(val_sal))];

    b = bar(1:2, m_bars, 'FaceColor','flat', 'EdgeColor','none', 'BarWidth',0.6);
    b.CData(1,:) = col_oxy;
    b.CData(2,:) = [0.7 0.7 0.7];

    errorbar(1:2, m_bars, s_bars, 'k.', 'LineWidth',1, 'CapSize',6);

    jitter = 0.15;
    scatter(ones(size(val_oxy)) + (rand(size(val_oxy))-0.5)*jitter, val_oxy, ...
        8, 'r', 'filled', 'MarkerFaceAlpha',0.3);
    scatter(2*ones(size(val_sal)) + (rand(size(val_sal))-0.5)*jitter, val_sal, ...
        8, 'k', 'filled', 'MarkerFaceAlpha',0.3);

    if numel(val_oxy) >= 3 && numel(val_sal) >= 3
        [~, p_bar] = ttest2(val_oxy, val_sal);
    else
        p_bar = NaN;
    end

    max_val = max([val_oxy; val_sal], [], 'omitnan');
    if isempty(max_val) || ~isfinite(max_val) || max_val <= 0
        max_val = 0.1;
    end

    y_sig_line = max_val * 1.15;
    line([1 2], [y_sig_line y_sig_line], 'Color','k', 'LineWidth',1);

    sig_txt = 'n.s.';
    if ~isnan(p_bar)
        if p_bar < 0.001
            sig_txt = '***';
        elseif p_bar < 0.01
            sig_txt = '**';
        elseif p_bar < 0.05
            sig_txt = '*';
        end
    end

    text(1.5, y_sig_line, sig_txt, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12, ...
        'FontWeight','normal');

    set(gca, 'Color','none', ...
        'Box','off', ...
        'XTick',[1 2], ...
        'XTickLabel',{'oxy','sal'}, ...
        'FontSize',8, ...
        'FontWeight','normal', ...
        'LineWidth',1, ...
        'TickDir','out');

    xlim([0.3 2.7]);
    % ylim([0, max_val*1.5]);
    title('5–20 min', 'FontSize',8, 'FontWeight','normal');
end

% sgtitle('Old-style time-resolved low-\theta analysis', 'FontWeight','bold');