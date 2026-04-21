%%
addpath("./main/stat")
% addpath utils
% load("./data/processed/psd_split_20230824.mat")
load("./data/processed/psd_split_20240416_1min.mat")
%load("./data_prev/skip_eeg_id.mat")
% 
% for i = 1:length(skip_eeg_id)
%     is_oxytocin(skip_eeg_id(i)) = -1;
% end

%% get p-value (Wilcoxon rank sum test)
zpsd_split_set = (psd_split_set - m_psd_set) ./ s_psd_set;
spec_sal = zpsd_split_set(:,:,:,is_oxytocin==0);
spec_oxy = zpsd_split_set(:,:,:,is_oxytocin==1);

sz = size(spec_sal);
pvals = nan(sz(1:3));

pbar = progressbar(sz(1), "stat test");
for nr = 1:sz(1)
    for nc = 1:sz(2)
        for nch = 1:3
            x = squeeze(spec_sal(nr, nc, nch, :));
            y = squeeze(spec_oxy(nr, nc, nch, :));

            % remove nan
            x = x(~isnan(x)); y = y(~isnan(y));
            if isempty(x) || isempty(y)
                continue;
            end

            pvals(nr, nc, nch) = ranksum(x, y);
        end
    end
    pbar.update(nr)
end
%
% %%
% alpha = 0.05;
%
% ch  = 3;
%
% figure;
% alpha = double(pvals(:,:,ch) < 0.05);
% alpha(alpha == 0) = 1;
%
% im = median(spec_oxy(:,:,ch,:), 4, "omitnan") - median(spec_sal(:,:,ch,:), 4, "omitnan");
% % im = median(spec_oxy(:,:,ch,:), 4, "omitnan");
%
% imagesc(tpsd/60, fpsd, im, "alphadata", alpha)
% colormap jet
% colorbar
% caxis([-.7, .7])
% ylim([1, 60])
% axis xy
% %
%% Plotting with alpha for significant regions
ch_names = {'mPFC', 'BLA', 'A1'};
alpha = 0.05; % significance level

fig = figure(1);
set(fig, 'Color', [1 1 1], "Position",  [620 748 1181 219])
for ch = 1:3 % choose channel

    % Calculate the median difference image
    im = median(spec_oxy(:,:,ch,:), 4, "omitnan") - median(spec_sal(:,:,ch,:), 4, "omitnan");

    % Create alpha mask based on p-values
    alpha_mask = double(pvals(:,:,ch) < alpha);
    alpha_mask(alpha_mask == 0) = 0; % set non-significant regions to be partially transparent
    alpha_mask(alpha_mask == 1) = 1;   % keep significant regions fully opaque

    % Apply Gaussian Smoothing and Interpolation
    im_smooth = imgaussfilt(im, 1);  % Apply Gaussian filter with sigma=1 for smoothing
    alpha_smooth = imgaussfilt(alpha_mask, 1);  % Apply Gaussian filter on alpha mask

    % Increase resolution for smoother display
    interp_factor = 3;  % Interpolation factor (adjust for more/less smoothness)
    [tpsd_interp, fpsd_interp] = meshgrid(linspace(min(tpsd/60), max(tpsd/60), size(im, 2) * interp_factor), ...
    linspace(min(fpsd), max(fpsd), size(im, 1) * interp_factor));

    im_interp = interp2(tpsd/60, fpsd, im_smooth, tpsd_interp, fpsd_interp, 'linear');
    alpha_interp = interp2(tpsd/60, fpsd, alpha_smooth, tpsd_interp, fpsd_interp, 'linear');

    % Plot the smoothed and interpolated data
    subplot(1,3,ch)
    i = imagesc(tpsd_interp(1,:), fpsd_interp(:,1), im_interp, 'AlphaData', alpha_interp);
set(gca, "LineWidth",1, ...
    'XTick', [-19.5 0 20 40 59.5], ...
    'XTickLabel', {'-20','0','20','40', '60'},...
    'YTick', [ 35 70]);
    
    colormap jet;
    cb = colorbar; ylabel(cb, 'z-scored power')
    cb.LineWidth = 1;
    clim([-.6, .6]);
    cb.Ticks = [-0.6 0 0.6];
    ylim([1, 70]);
    title(sprintf('%s', ch_names{ch}))
    xlabel('Time (min)');
    ylabel('Frequency (Hz)');
   % title(sprintf('Smoothed Significant Regions (p < 0.05) - Channel %d', ch));

    % Add vertical dotted line at x = 0
    xline(0, '--k', 'LineWidth', 1); % Dotted black line at x=0

    % Additional plot adjustments
    axis xy;


end
