clear; clc; close all;


data_file   = './data/ERP_all_temporal.mat';
meta_file   = './data/meta_20251011.mat';
info_file   = './data/final_0428a-trialInfo.mat';

exclude_sessions = [25, 39, 55, 56, 57, 58, 59, 60, 67, 68, 69, 70];
conds_to_analyze = [1 2 3 4];   
n_mean      = 10;          % 블록당 트라이얼 수
srate_onset = 300000;      


fprintf('Loading data...\n');
load(data_file, 'ERP_all', 'drug', 'times', 'srate', 'chan_names');
load(meta_file, 'meta');
load(info_file, 'trialInfo');
chan_names = {'mPFC', 'BLA', 'AC'};

trial_list_universal = single(trialInfo.trial_list(:));
trial_onset_universal = single(trialInfo.trial_onset(:));
pup_indices = find(ismember(trial_list_universal', conds_to_analyze));


n_total = size(ERP_all, 4);
valid_mask = true(1, n_total);
valid_mask(exclude_sessions) = false;

groups(1).name = 'oxytocin'; 
groups(1).mask = (drug == 1) & valid_mask;
groups(1).ids  = find(groups(1).mask);
groups(2).name = 'saline';   
groups(2).mask = (drug ~= 1) & valid_mask;
groups(2).ids  = find(groups(2).mask);

groups(1).color = [0.85 0.20 0.2]; % Red
groups(2).color = [0.3 0.3 0.3];   % Gray

N_match = min(length(groups(1).ids), length(groups(2).ids));
onset_sec_all = trial_onset_universal(pup_indices) / srate_onset;

%% Supple Figure 1A: Stacked ERP in mPFC, BLA
chs_to_plot = [1 1 2 2];   % mPFC, mPFC, BLA, BLA
conds = [1 2 1 2];         % oxy, sal, oxy, sal
positions = {
    [0.1300 0.5461 0.1466 0.4212],  % mPFC oxy
    [0.3075 0.5461 0.1466 0.4212],  % mPFC sal
    [0.4850 0.5461 0.1466 0.4212],  % BLA oxy
    [0.6625 0.5461 0.1466 0.4212]   % BLA sal
    };

ax_all = gobjects(1,4);

for i = 1:4
    ch = chs_to_plot(i);
    cond = conds(i);

    ax = subplot(5,4,i); % dummy (Position으로 덮어씀)
    ax.Position = positions{i};
    hold(ax,'on');

    for blk = 1:Ng
        offset = (blk-1)*offset_val;
        plot(ax, times, groups(cond).data(ch,:,blk)-offset, ...
            'Color', cmap(min(blk,Ng),:), 'LineWidth', 1);
    end

    xline(ax,0,'--k','LineWidth',1);

    title(ax, sprintf('%s | %s', chan_names{ch}, groups(cond).name));

    xlim(ax, plot_time_range);
    ylim(ax, [-670 50]);

    box(ax,'off');
    set(ax, "XTick", xticks, 'YTick', yticks, ...
        'YTickLabel', yticklabels, 'YDir','normal', ...
        'LineWidth',1, 'FontSize',10);

    set(ax,'XColor','none','YColor','none');

    ax_all(i) = ax;
end

%% 8. Peak Amplitude Tracking (Mean of Window)
fprintf('Tracking Window-Mean Amplitude over time...\n');


peak_windows = {[0.005 0.035], [0.005 0.035], [0.005 0.035]}; 
stat_alpha   = 0.05;

fig_trend = figure(20); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1100 300]);
plot_chs = [3 1 2]; 

for i = 1:3
    ch_idx = plot_chs(i);
    win = peak_windows{ch_idx};
    win_idx = (times >= win(1)) & (times <= win(2));
    
    subplot(1, 3, i); hold on;
    
    group_data = cell(1, 2);
    block_times_min = [];
    h_line = gobjects(1,2);

    for g = 1:2
        ids = groups(g).ids(1:N_match);
        n_subj = length(ids);
        
        % 세션별 블록 데이터 계산
        for s = 1:n_subj
            subj_dat = squeeze(ERP_all(ch_idx, :, pup_indices, ids(s)));
            nTrials = size(subj_dat, 2);
            n_blocks = floor(nTrials / n_mean);

            if s == 1
                subj_block_peaks = nan(n_subj, n_blocks);
                if isempty(block_times_min), block_times_min = nan(1, n_blocks); end
            end

            for b = 1:n_blocks
                idx = (b-1)*n_mean + (1:n_mean);
                trace = mean(subj_dat(:, idx), 2, 'omitnan');
                roi_data = trace(win_idx);
                
                % [수정] 노이즈에 강건하도록 Peak가 아닌 윈도우 내 '평균(Mean)' 사용
                subj_block_peaks(s, b) = mean(roi_data, 'omitnan');

                if g == 1 && s == 1
                    block_times_min(b) = mean(onset_sec_all(idx), 'omitnan') / 60;
                end
            end
        end
        group_data{g} = subj_block_peaks;


        m_vals = mean(subj_block_peaks, 1, 'omitnan');
        s_vals = std(subj_block_peaks, 0, 1, 'omitnan') ./ sqrt(n_subj);
        lineProps = {'Color', groups(g).color, 'LineWidth', 1.5};

        if exist('shadedErrorBar', 'file') == 2
            h_sh = shadedErrorBar(block_times_min, m_vals, s_vals, lineProps, 1);
            h_line(g) = h_sh.mainLine;


            if isfield(h_sh, 'patch') && all(isgraphics(h_sh.patch))
                set(h_sh.patch, 'EdgeColor', 'none', 'FaceAlpha', 0.2); 
            end
            if isfield(h_sh, 'edge') && all(isgraphics(h_sh.edge))
                set(h_sh.edge, 'Visible', 'off'); 
            end
        else
            h_line(g) = plot(block_times_min, m_vals, 'Color', groups(g).color, 'LineWidth', 1.5);
        end
    end


    n_blocks = length(block_times_min);
    sig_blocks = false(1, n_blocks);
    for b = 1:n_blocks
        v1 = group_data{1}(:, b); v2 = group_data{2}(:, b);
        v1_c = v1(~isnan(v1)); v2_c = v2(~isnan(v2));
        if length(v1_c) >= 3 && length(v2_c) >= 3
            sig_blocks(b) = (ranksum(v1_c, v2_c) < stat_alpha);
        end
    end

    
    title(sprintf('%s (%.0f-%.0f ms Mean)', chan_names{ch_idx}, win(1)*1000, win(2)*1000));
    xlabel('Time from injection (min)');
    if i == 1, ylabel('Mean Amplitude (\muV)'); end
    
    set(gca, 'LineWidth', 1, 'FontSize', 11, 'Box', 'off', 'TickDir', 'in', 'YDir', 'reverse');
    grid off; axis square;

    if i == 1, ylim([-140 20]); else, ylim([-80 20]); end
    
    
    y_limits = get(gca, 'YLim');
    line_y = y_limits(1) * 0.92; 
    sig_idx = find(sig_blocks);
    if ~isempty(sig_idx)
        d = diff(sig_idx);
        split_pts = [0 find(d > 1) length(sig_idx)];
        for k = 1:length(split_pts)-1
            idx_range = sig_idx(split_pts(k)+1 : split_pts(k+1));
            plot([block_times_min(idx_range(1)) block_times_min(idx_range(end))], ...
                 [line_y line_y], 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
    end
    
    if i == 1
        lgd = legend(h_line, {groups(1).name, groups(2).name}, 'Location', 'southwest', 'Box', 'off');
        lgd.ItemTokenSize = [5, 12]; % 사용자님 선호 스타일
    end
end
fprintf('Tracking Analysis Done.\n');
