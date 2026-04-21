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


%% 실행 함수 호출
gc_lineplot_with_inset(gc_all, drug);

function gc_lineplot_with_inset(gc_all, drug)
% --- Parameter Setup ---
freq = 1:55;
theta_range = [8 12];

chanLabel = {'mPFC', 'BLA', 'AC'};
chanComb = [1 2; 1 3; 2 3];
timeIdx = 2;

% --- Outlier Removal ---
sal_outliers = [4, 13, 14, 16, 18, 19];
oxy_outliers = [2, 16];

sal_idx = find(drug == 0); sal_idx(sal_outliers) = [];
oxy_idx = find(drug == 1); oxy_idx(oxy_outliers) = [];

% --- Main Loop ---
for combIdx = 2 % mPFC-A1 조합
    for dirIdx = 1:2
        % ---------------------------------------------------
        % [A] Main Subplot (Line Plot)
        % ---------------------------------------------------
        % sub_ax = subplot(2, 4, spIdx); 
        sub_ax = subplot(2,3,-dirIdx+4); hold on; axis square
        cla; hold on;

        % Data Extraction
        gc_sal_raw = squeeze(gc_all(:, combIdx, dirIdx, sal_idx, timeIdx))';
        gc_oxy_raw = squeeze(gc_all(:, combIdx, dirIdx, oxy_idx, timeIdx))';

        valid_sal = all(~isnan(gc_sal_raw), 2) & any(gc_sal_raw ~= 0, 2);
        valid_oxy = all(~isnan(gc_oxy_raw), 2) & any(gc_oxy_raw ~= 0, 2);

        gc_sal = gc_sal_raw(valid_sal, freq);
        gc_oxy = gc_oxy_raw(valid_oxy, freq);

        mean_sal = nanmedian(gc_sal, 1);
        sem_sal = nanstd(gc_sal, 0, 1) / sqrt(size(gc_sal,1));
        mean_oxy = nanmedian(gc_oxy, 1);
        sem_oxy = nanstd(gc_oxy, 0, 1) / sqrt(size(gc_oxy,1));

        ms_sal = smoothdata(mean_sal, 'movmean', 9);
        ss_sal = smoothdata(sem_sal, 'movmean', 5);
        ms_oxy = smoothdata(mean_oxy, 'movmean', 9);
        ss_oxy = smoothdata(sem_oxy, 'movmean', 5);

        % Main Plotting
        fill([freq, fliplr(freq)], [ms_sal+ss_sal, fliplr(ms_sal-ss_sal)], ...
            [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(freq, ms_sal, 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

        fill([freq, fliplr(freq)], [ms_oxy+ss_oxy, fliplr(ms_oxy-ss_oxy)], ...
            [1 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(freq, ms_oxy, 'Color', [0.85 0.2 0.2], 'LineWidth', 1);

        % Main Stars
        for f = 1:length(freq)
            [~, p] = ttest2(gc_sal(:,f), gc_oxy(:,f));
            if p < 0.05
                col = 'k';
                if ms_oxy(f) > ms_sal(f), col = 'r'; end
                text(freq(f), 0.23, '*', 'HorizontalAlignment', 'center', ...
                    'FontSize', 10, 'Color', col, 'FontWeight', 'normal');
            end
        end

        % --- 축 설정 수정 (Tick 방향 및 레이블 제거) ---
        xlim([1 55]);
        xticks([10 20 30 40 50]);
        ylim([0.05 0.25]);
        yticks([0.05 0.10 0.15 0.20 0.25]);

        % [수정] TickDir을 'in'으로 설정
        set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 1, ...
            'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'normal');

        % [수정] 두 번째 subplot(오른쪽)의 Y축 숫자 제거
        if dirIdx == 2
            % set(gca, 'YTickLabel', []);
        end

        from_ch = chanLabel{chanComb(combIdx, dirIdx)};
        to_ch   = chanLabel{chanComb(combIdx, 3 - dirIdx)};
        title(sprintf('%s \\rightarrow %s', from_ch, to_ch), ...
            'FontWeight', 'normal', 'FontSize', 10);

        % [수정] ylabel 텍스트를 Frequency (Hz)로 변경

            ylabel('GC (a.u.)', 'FontWeight', 'normal');

            xlabel('Freq (Hz)', 'FontWeight', 'normal');

        % ---------------------------------------------------
        % [B] Inset Bar Plot (기존 동일)
        % ---------------------------------------------------
        pos = get(sub_ax, 'Position');
        inset_pos = [pos(1)+pos(3)*0.55, pos(2)+pos(4)*0.35, pos(3)*0.35, pos(4)*0.35];

        axes('Position', inset_pos); hold on;
        % ... (이하 Inset 코드 생략) ...

        f_idx = freq >= theta_range(1) & freq <= theta_range(2);
        val_sal = mean(gc_sal(:, f_idx), 2);
        val_oxy = mean(gc_oxy(:, f_idx), 2);

        m_bars = [mean(val_oxy), mean(val_sal)];
        s_bars = [std(val_oxy)/sqrt(length(val_oxy)), std(val_sal)/sqrt(length(val_sal))];

        b = bar(1:2, m_bars, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.6);
        b.CData(2,:) = [0.7 0.7 0.7];
        b.CData(1,:) = [0.85 0.2 0.2];

        errorbar(1:2, m_bars, s_bars, 'k.', 'LineWidth', 1, 'CapSize', 6);

        jitter = 0.15;
        scatter(ones(size(val_oxy)) + (rand(size(val_oxy))-0.5)*jitter, val_oxy, ...
            8, 'r', 'filled', 'MarkerFaceAlpha', 0.3);
        scatter(2*ones(size(val_sal)) + (rand(size(val_sal))-0.5)*jitter, val_sal, ...
            8, 'k', 'filled', 'MarkerFaceAlpha', 0.3);

        % --- Inset Statistics ---
        [~, p_bar] = ttest2(val_oxy, val_sal);

        max_val = max([val_oxy; val_sal]);
        if max_val <= 0, max_val = 0.1; end

        y_sig_line = max_val * 1.15;
        line([1 2], [y_sig_line y_sig_line], 'Color', 'k', 'LineWidth', 1);

        sig_txt = 'n.s.';
        if p_bar < 0.001, sig_txt = '***';
        elseif p_bar < 0.01, sig_txt = '**';
        elseif p_bar < 0.05, sig_txt = '*';
        end

        if p_bar < 0.05
            text(1.5, y_sig_line, sig_txt, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'normal');
        else
            text(1.5, y_sig_line, sig_txt, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'normal');
        end

        % --- Inset Styling ---
        set(gca, 'Color', 'none', 'Box', 'off', 'XTick', [1 2], ...
            'XTickLabel', {'oxy', 'sal'}, 'FontSize', 10, 'FontWeight', 'normal', ...
            'LineWidth', 1); % LineWidth 1 고정

        xlim([0.3 2.7]);
        ylim([0, max_val * 1.5]); % Inset Y축 넉넉하게

        title('\theta_{low}', 'FontSize', 10, 'FontWeight', 'normal');
    end
end
end






%% driving synchrony, gating information
% 1. 초기 설정 및 파라미터
% clearvars -except psd gc plv drug meta inj_time; % 데이터 로드 상태 유지

% 데이터 로드 (필요시)
if ~exist('psd', 'var'), load('data/processed/psd_time_250511d.mat'); end
if ~exist('gc', 'var'),  load('data/processed/gc_cleaned.mat'); end
if ~exist('plv', 'var'), load('./data/processed/plv_0401b.mat'); end
if ~exist('drug_type', 'var'), load('data/drug.mat'); end
if ~exist('meta', 'var'), meta = readtable('./data/meta.xlsx'); end
load('/Users/dayoung/Documents/Dayoung/Dayoung2023/resting/data/processed/plv_0401b.mat', 'inj_time')



%% 1. 파라미터 및 데이터 추출 루프
% 채널 설정
ch_mPFC = 1;  % mPFC 채널
% ch_BLA   = 2;  % BLA
ch_A1    = 3;  % AC (A1) 채널
dir = 2;

% 분석 파라미터
n_sessions = 124;
duration_post = 15 * 60; % 주입 후 20분
range_psd = [8 12];      % Beta/Gamma range
range_plv = [8 12];       % Theta low
range_gc  = [8 12];       % Theta low

fprintf('데이터 추출 중 (mPFC 채널 오류 수정 완료)... \n');
data_mat = nan(n_sessions, 5); % [AC_P, mPFC_P, PLV, GC, Drug]

for i = 1:n_sessions
    % 유효성 체크 (데이터 존재 여부 및 Injection Time 확인)
    if i > size(psd,1) || isempty(psd{i,1}) || isempty(plv{i,1}) || isempty(gc{i,1}) || isnan(inj_time(i,2))
        continue;
    end
    
    % --- 시간 윈도우 설정 (주입 1분 후 ~ 21분) ---
    t_start = inj_time(i, 2) + 5*60; 
    t_end   = t_start + duration_post;
    
    % 각 데이터별 타임 인덱스 추출
    idx_p = (psd{i,2} > t_start) & (psd{i,2} <= t_end);
    idx_l = (plv{i,2} > t_start) & (plv{i,2} <= t_end);
    idx_g = (gc{i,2} > t_start) & (gc{i,2} <= t_end);
    
    % if ~any(idx_p) || ~any(idx_l) || ~any(idx_g), continue; end
    
    % --- 1. AC(A1) Power 계산 ---
    f_p = linspace(0, 121, size(psd{i,1}, 2));
    c_p = (f_p >= range_psd(1)) & (f_p <= range_psd(2));
    raw_p_ac = psd{i,1}(idx_p, c_p, ch_A1);
    val_psd_AC = 10 * log10(mean(raw_p_ac(:), 'omitnan'));
    
    % --- 2. mPFC Power 계산 (채널 인덱스 수정) ---
    % 기존 코드에서 ch_A1을 쓰고 있던 부분을 ch_mPFC로 수정했습니다.
    raw_p_mpfc = psd{i,1}(idx_p, c_p, ch_mPFC);
    val_psd_mPFC = 10 * log10(mean(raw_p_mpfc(:), 'omitnan'));
    
    % --- 3. PLV (AC-mPFC Sync) ---
    f_l = linspace(0, 100, size(plv{i,1}, 2));
    c_l = (f_l >= range_plv(1)) & (f_l <= range_plv(2));
    raw_l = plv{i,1}(idx_l, c_l, 2); % ch_plv=2 (AC-mPFC 조합)
    val_plv = mean(raw_l(:), 'omitnan');
    
    % --- 4. GC (AC -> mPFC Gating) ---
    f_g = linspace(0, 512, size(gc{i,1}, 2));
    c_g = (f_g >= range_gc(1)) & (f_g <= range_gc(2));
    raw_g = gc{i,1}(idx_g, c_g, 2, dir); % AC -> mPFC 방향
    val_gc = mean(raw_g(:), 'omitnan');
    
    % 데이터 매트릭스 저장
    data_mat(i, :) = [val_psd_AC, val_psd_mPFC, val_plv, val_gc, drug(i)];
end

% 2. 데이터 클리닝 및 필터링
% 1) 유효하지 않은 값(NaN, Inf) 제거
idx_invalid = any(~isfinite(data_mat), 2);

% 2) 비정상적인 GC 값 필터링 (0 근처 혹은 0.3 초과)
idx_gc_outlier = (data_mat(:, 4) <= 1e-5) | (data_mat(:, 4) > 0.3);

% 3) Power Outlier (실험 환경 노이즈로 인한 비정상 수치 제거)
% mPFC나 AC 파워가 지나치게 낮은 세션 제외
idx_low_power = (data_mat(:, 1) < -20) | (data_mat(:, 2) < -20); 

% 최종 마스크 생성
valid_mask = ~idx_invalid & ~idx_gc_outlier & ~idx_low_power;
clean_data = data_mat(valid_mask, :);

fprintf('--- 필터링 요약 ---\n');
fprintf('전체 세션: %d -> 최종 분석 대상: %d\n', n_sessions, size(clean_data, 1));
fprintf('제외 사유: 유효성(%d), GC이상치(%d), Power저하(%d)\n', ...
    sum(idx_invalid), sum(idx_gc_outlier & ~idx_invalid), sum(idx_low_power & ~idx_invalid & ~idx_gc_outlier));
% 상관관계 분석 플롯 (Panel 5, 6, 7, 8) - 텍스트 오른쪽 정렬 및 스타일 최종본
plot_configs = {
    'Driving Sync (AC)',   1, 3, 'AC Power (dB)',   'AC-mPFC \theta_{high} PLV', 1;
    'Driving Sync (mPFC)', 2, 3, 'mPFC Power (dB)', 'AC-mPFC \theta_{high} PLV', 2;
    'Gating Info (AC)',    1, 4, 'AC Power (dB)',   'mPFC\rightarrowAC \theta_{high} GC',  3;
    'Gating Info (mPFC)',  2, 4, 'mPFC Power (dB)', 'mPFC\rightarrowAC \theta_{high} GC',  4
};

g_oxy = clean_data(clean_data(:,5) == 1, :); 
g_sal = clean_data(clean_data(:,5) == 0, :); 
colors = {[0.85 0.2 0.2], [0.3 0.3 0.3]}; 
alphas = [0.4, 0.3];
std_font_sz = 10; 
fig = figure(1);clf;
set(fig, "Position", [700 965 1040 273], "Color" , 'w')

for i = 1:4
    target_pos = plot_configs{i, 6};
    subplot(1, 4, target_pos); 
    cla; % 이전 설정 초기화
    axis square; hold on;
    
    title_str = plot_configs{i, 1};
    x_idx = plot_configs{i, 2};
    y_idx = plot_configs{i, 3};
    
    groups = {g_oxy, g_sal};
    r_vals = [NaN, NaN]; p_vals = [NaN, NaN];
    h_leg_proxy = []; 
    
    % 1. X축 범위 고정
    xlim([2 10]);
    ax_limit = xlim;

    % 2. 스캐터 플롯
    for g = 1:2
        curr_data = groups{g};
        scatter(curr_data(:, x_idx), curr_data(:, y_idx), 20, colors{g}, ...
            'filled', 'MarkerFaceAlpha', alphas(g), 'HandleVisibility', 'off');
    end
    
    % 3. 회귀 직선 및 통계
    for g = 1:2
        curr_data = groups{g};
        x_val = curr_data(:, x_idx);
        y_val = curr_data(:, y_idx);
        
        if length(x_val) > 2
            [r, p] = corr(x_val, y_val, 'Rows', 'complete');
            r_vals(g) = r; p_vals(g) = p;
            
            l_style = ':'; l_width = 1.2;
            if p < 0.05, l_style = '-'; l_width = 1.8; end
            
            p_fit = polyfit(x_val, y_val, 1);
            x_range = linspace(ax_limit(1), ax_limit(2), 100); 
            plot(x_range, polyval(p_fit, x_range), 'Color', colors{g}, ...
                'LineStyle', l_style, 'LineWidth', l_width, 'HandleVisibility', 'off');
            
            % 범례용 프록시 (항상 실선)
            h_leg_proxy(g) = plot(nan, nan, 'Color', colors{g}, 'LineStyle', '-', 'LineWidth', 1.5);
        end
    end
    
    % --- 축 스타일 설정 (Tick 'in', Y-label 복구) ---
    set(gca, 'Box', 'off', 'TickDir', 'in', 'LineWidth', 1, 'FontSize', std_font_sz);
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto'); % 레이블 강제 활성화
    
    xlabel(plot_configs{i, 4}, 'FontSize', std_font_sz);
    ylabel(plot_configs{i, 5}, 'FontSize', std_font_sz);
    % title(title_str, 'FontSize', std_font_sz + 1, 'FontWeight', 'normal');
    
    % --- 통계 텍스트 (오른쪽 정렬 수정) ---
    text(0.95, 0.94, sprintf('oxy: r=%.2f, p=%.2f', r_vals(1), p_vals(1)), ...
        'Units', 'normalized', 'Color', colors{1}, 'FontSize', 10, ...
        'HorizontalAlignment', 'right', 'FontWeight', 'normal');
    text(0.95, 0.86, sprintf('sal: r=%.2f, p=%.2f', r_vals(2), p_vals(2)), ...
        'Units', 'normalized', 'Color', colors{2}, 'FontSize', 10, ...
        'HorizontalAlignment', 'right', 'FontWeight', 'normal');

    % --- 범례 설정 (토큰 짧게) ---
    % if target_pos == 8
    %     [lgd, ~] = legend(h_leg_proxy, {'oxytocin', 'saline'}, 'Location', 'southeast', 'Box', 'off');
    %     lgd.FontSize = 8;
    %     lgd.ItemTokenSize = [10, 8]; 
    % end
end

%%
fig = fig500;
exportgraphics(fig, 'Main_Figure3_260203a.pdf', 'ContentType', 'vector');
exportgraphics(fig, 'Main_Figure3_260203a.tif', 'Resolution', 600);
