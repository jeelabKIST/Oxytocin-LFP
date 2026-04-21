%% GC plot
% addpath hb_sptools/
load('data/processed/gc_cleaned.mat')

inj_time = meta_info.inj_time;
drug = meta_info.drug;
skip_ratio = meta_info.skip_ratio;

% 유효한 파일 인덱스 찾기
isEmptyGC = cellfun(@isempty, gc(:,1));
valid_fileIdx = find((skip_ratio <= 10) & ~isEmptyGC);
drug = drug(valid_fileIdx);
inj_time = inj_time(valid_fileIdx,:);

% GC 데이터 초기화
gc_all = nan([512, 3, 2, length(valid_fileIdx), 3]);
for newIdx = 1:length(valid_fileIdx)
    fileIdx = valid_fileIdx(newIdx);
    gc_data = real(gc{fileIdx,1});
    t = gc{fileIdx,2};

    if inj_time(newIdx,1) > 0
        inj_start_end = hb_findIdx(inj_time(newIdx,:), t);
        for timeIdx = 1:3
            switch(timeIdx)
                case 1
                    time_win = inj_start_end(1) - 20*60 + 1 : inj_start_end(1);
                    time_win = time_win(time_win > 0);
                case 2
                    time_win = inj_start_end(end) : inj_start_end(end) + 20*60;
                case 3
                    time_win = inj_start_end(end) + 20*60 : length(t);
            end
            indiv_gc = nanmedian(gc_data(time_win,:,:,:), 1);
            gc_all(:,:,:,newIdx,timeIdx) = indiv_gc;
        end
    end
end

oxy  = find(drug == 1);
sal  = find(drug == 0);




%% [1] 통합 피규어 초기 설정
fig = figure(2); clf;
set(fig, 'Color', 'w', 'Position', [269 517 918 506]);

% 공통 스타일 파라미터
c_oxy = [0.85 0.2 0.2]; c_sal = [0.3 0.3 0.3];
std_lw = 1; std_fs = 10;
freq_line = 1:55; 
freq_bar = 1:121;
freq_bands = {[4 8], [8 12], [18 24], [24 32], [35 50], [70 90]};
band_names = {'\theta_{low}', '\theta_{high}', '\beta_{low}', '\beta_{high}', '\gamma_{low}', '\gamma_{high}'};
chan_labels = {'mPFC', 'BLA', 'AC'};
chan_comb = [1 2; 1 3; 2 3]; % 1:mPFC-BLA, 2:mPFC-AC, 3:BLA-AC

% 아웃라이어 제거된 인덱스 (기존 logic 유지)
sal_idx = find(drug == 0); sal_idx([4, 13, 14, 16, 18, 19]) = [];
oxy_idx = find(drug == 1); oxy_idx([2, 16]) = [];


% [2] Main Plotting Loop (3 Rows x 4 Columns)
for combIdx = 1:3 
    for dirIdx = 1:2 
        
        % 데이터 추출 (timeIdx=2: Post-injection)
        raw_sal = squeeze(gc_all(:, combIdx, dirIdx, sal_idx, 2))';
        raw_oxy = squeeze(gc_all(:, combIdx, dirIdx, oxy_idx, 2))';
        
        % ---------------------------------------------------
        % [A] Line Plot (인덱스: 1, 2 / 5, 6 / 9, 10)
        % ---------------------------------------------------
        sp_line = (combIdx-1)*4 + dirIdx; 
        subplot(3, 4, sp_line); hold on; axis square;
        
        % 통계 및 스무딩
        ms_sal = smoothdata(nanmean(raw_sal(:, freq_line), 1), 'movmean', 9);
        ss_sal = smoothdata(nanstd(raw_sal(:, freq_line), 0, 1)/sqrt(size(raw_sal,1)), 'movmean', 5);
        ms_oxy = smoothdata(nanmean(raw_oxy(:, freq_line), 1), 'movmean', 9);
        ss_oxy = smoothdata(nanstd(raw_oxy(:, freq_line), 0, 1)/sqrt(size(raw_oxy,1)), 'movmean', 5);
        
        % 그림 그리기 (Shaded Area & Line)
        fill([freq_line, fliplr(freq_line)], [ms_sal+ss_sal, fliplr(ms_sal-ss_sal)], c_sal, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(freq_line, ms_sal, 'Color', c_sal, 'LineWidth', std_lw);
        fill([freq_line, fliplr(freq_line)], [ms_oxy+ss_oxy, fliplr(ms_oxy-ss_oxy)], c_oxy, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(freq_line, ms_oxy, 'Color', c_oxy, 'LineWidth', std_lw);
        
        % --- [추가] Line Plot 주파수별 통계 (ttest2) ---
        for f = 1:length(freq_line)
            [~, p_val] = ttest2(raw_sal(:, f), raw_oxy(:, f));
            if p_val < 0.05
                star_col = 'k'; % 기본 검정색
                if ms_oxy(f) > ms_sal(f), star_col = c_oxy; end % Oxy > Sal일 때 빨간색
                text(freq_line(f), 0.23, '*', 'HorizontalAlignment', 'center', ...
                    'FontSize', 10, 'Color', star_col, 'FontWeight', 'normal');
            end
        end
        
        % 스타일 설정
        set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', std_lw, 'FontSize', std_fs, 'XLim', [1 55], 'YLim', [0.05 0.25]);
        from_node = chan_labels{chan_comb(combIdx, dirIdx)};
        to_node = chan_labels{chan_comb(combIdx, 3-dirIdx)};
        title(sprintf('%s \\rightarrow %s', from_node, to_node), 'FontSize', std_fs, 'FontWeight', 'normal');
        if dirIdx == 1, ylabel('GC (a.u.)'); end
        if combIdx == 3, xlabel('Freq (Hz)'); end

        % ---------------------------------------------------
        % [B] Bar Graph (인덱스: 3, 4 / 7, 8 / 11, 12)
        % ---------------------------------------------------
        sp_bar = (combIdx-1)*4 + dirIdx + 2;
        subplot(3, 4, sp_bar); hold on;
        
        m_vals = []; s_vals = []; p_vals = [];
        for b = 1:6
            f_idx = hb_findIdx(freq_bands{b}, freq_bar);
            val_sal = nanmean(raw_sal(:, f_idx), 2);
            val_oxy = nanmean(raw_oxy(:, f_idx), 2);
            m_vals(b, :) = [nanmean(val_sal), nanmean(val_oxy)];
            s_vals(b, :) = [nanstd(val_sal)/sqrt(sum(~isnan(val_sal))), nanstd(val_oxy)/sqrt(sum(~isnan(val_oxy)))];
            p_vals(b) = ranksum(val_sal, val_oxy); 
        end
        
        % Bar Plotting
        x = 1:6; bw = 0.35;
        b1 = bar(x - bw/2, m_vals(:,1), bw, 'FaceColor', c_sal, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        b2 = bar(x + bw/2, m_vals(:,2), bw, 'FaceColor', c_oxy, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        errorbar(x - bw/2, m_vals(:,1), s_vals(:,1), 'k.', 'LineWidth', 0.8, 'CapSize', 0);
        errorbar(x + bw/2, m_vals(:,2), s_vals(:,2), 'k.', 'LineWidth', 0.8, 'CapSize', 0);
        
        % --- [수정] Bar Graph 유의성 표시 (조건부 색상) ---
        for b = 1:6
            if p_vals(b) < 0.05
                star_col = 'k';
                if m_vals(b, 2) > m_vals(b, 1), star_col = c_oxy; end % Oxy > Sal일 때 빨간색
                text(b, max(m_vals(b,:)) + s_vals(b, 1) + 0.02, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', star_col);
            end
        end
        
        % 스타일 설정
        set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', std_lw, 'FontSize', std_fs, 'XTick', x, 'XTickLabel', band_names, 'YLim', [0 0.25]);
        xtickangle(45);
        title(sprintf('%s \\rightarrow %s', from_node, to_node), 'FontSize', std_fs, 'FontWeight', 'normal');
        if dirIdx == 1, ylabel('GC (a.u.)'); end
    end
end

%%
exportgraphics(fig, 'Supple_Figure3_260130a.pdf', 'ContentType', 'vector');
exportgraphics(fig, 'Supple_Figure3_260130a.tif', 'Resolution', 600);


