function [act_matrix, exg_matrix, time_vec, fig] = plot_hdsemg_time_heatmap(emg, fs_emg, window_ms, title_str, overlap_ratio, normalize_baseline, baseline_s)
% PLOT_HDSEMG_TIME_HEATMAP  Heatmap "канал x время": RMS-амплитуда в
%   скользящих окнах длиной window_ms (мс) с перекрытием. Строит ДВА
%   subplot в одном окне: (1) HD-sEMG 126 каналов, (2) монополярные
%   EXG1-8. Подписывает по оси Y первый и последний электрод каждого
%   столбца сетки (A1/A9, B1/B9, ...) для HD-подграфика.
%
%   [act_matrix, exg_matrix, time_vec, fig] = PLOT_HDSEMG_TIME_HEATMAP(emg, fs_emg, window_ms, title_str, overlap_ratio, normalize_baseline, baseline_s)
%
%   Вход:
%     emg               - 134 x N массив EMG-сигналов (126 HD + 8 EXG)
%     fs_emg             - частота дискретизации EMG, Гц
%     window_ms          - длина окна в мс (по заданию = 128)
%     title_str           - заголовок (например, название движения)
%     overlap_ratio       - доля перекрытия окон, 0..1 (по умолчанию 0.5)
%     normalize_baseline  - true/false: нормировать каждый канал на его
%                           собственный уровень покоя (устраняет разницу
%                           в усилении/шумовом поле между каналами, видную
%                           как "радужная" полоса в состоянии покоя).
%                           По умолчанию true.
%     baseline_s          - длительность периода покоя в начале записи
%                           для расчёта базового уровня, с (по умолчанию 1.0)
%
%   Выход:
%     act_matrix - 126 x n_windows, RMS по HD-каналам (норм. или нет)
%     exg_matrix - 8 x n_windows, RMS по EXG-каналам (норм. или нет)
%     time_vec   - 1 x n_windows, время центра каждого окна, с
%     fig        - handle фигуры

    if nargin < 3 || isempty(window_ms)
        window_ms = 128;
    end
    if nargin < 5 || isempty(overlap_ratio)
        overlap_ratio = 0.5;
    end
    if nargin < 6 || isempty(normalize_baseline)
        normalize_baseline = false;   % по умолчанию ВЫКЛ - сырые RMS с обрезкой по перцентилям читаемее
    end
    if nargin < 7 || isempty(baseline_s)
        baseline_s = 1.0;
    end

    n_hd  = 126;
    n_exg = 8;
    emg_hd  = emg(1:n_hd, :);
    emg_exg = emg(n_hd+1:n_hd+n_exg, :);

    win_samples  = round(window_ms / 1000 * fs_emg);
    step_samples = max(1, round(win_samples * (1 - overlap_ratio)));
    n_samples    = size(emg, 2);
    starts       = 1:step_samples:(n_samples - win_samples + 1);
    n_windows    = numel(starts);

    act_matrix = zeros(n_hd, n_windows);
    exg_matrix = zeros(n_exg, n_windows);
    for w = 1:n_windows
        idx = starts(w) : starts(w) + win_samples - 1;
        act_matrix(:, w) = sqrt(mean(emg_hd(:, idx).^2, 2));
        exg_matrix(:, w) = sqrt(mean(emg_exg(:, idx).^2, 2));
    end
    time_vec = (starts + win_samples/2 - 1) / fs_emg;   % центр каждого окна, с

    % --- Нормализация по базовому уровню покоя (устраняет разницу в усилении между каналами) ---
    % Используем z-score (не отношение!): (x - baseline_mean) / baseline_std.
    % Отношение (x / baseline_mean) взрывается на тихих каналах с маленьким
    % baseline - крошечный шум делится на маленькое число и выглядит как
    % большая "активность". Z-score устойчив к этому.
    colorbar_label = 'RMS амплитуда';
    if normalize_baseline
        n_baseline_win = max(1, sum(time_vec <= baseline_s));

        baseline_hd_mean  = mean(act_matrix(:, 1:n_baseline_win), 2);
        baseline_hd_std   = std(act_matrix(:, 1:n_baseline_win), 0, 2);
        baseline_exg_mean = mean(exg_matrix(:, 1:n_baseline_win), 2);
        baseline_exg_std  = std(exg_matrix(:, 1:n_baseline_win), 0, 2);

        % защита от деления на ~0 (канал без шума в покое -> используем
        % медианный std по всем каналам как разумный минимум)
        min_std_hd  = max(median(baseline_hd_std),  eps);
        min_std_exg = max(median(baseline_exg_std), eps);
        baseline_hd_std(baseline_hd_std   < 0.1 * min_std_hd)  = min_std_hd;
        baseline_exg_std(baseline_exg_std < 0.1 * min_std_exg) = min_std_exg;

        act_matrix = (act_matrix - baseline_hd_mean) ./ baseline_hd_std;
        exg_matrix = (exg_matrix - baseline_exg_mean) ./ baseline_exg_std;
        colorbar_label = 'Отклонение от покоя, SD (z-score)';

        fprintf('  [baseline] z-score по %.2f с покоя (%d окон)\n', baseline_s, n_baseline_win);
    end

    % --- Подписи первого/последнего электрода каждого столбца (A1/A9...N1/N9) ---
    ytick_pos = zeros(1, 28);
    ytick_lab = cell(1, 28);
    k = 0;
    for c = 1:14
        col_letter = char('A' + c - 1);
        k = k + 1; ytick_pos(k) = (c-1)*9 + 1; ytick_lab{k} = sprintf('%s1', col_letter);
        k = k + 1; ytick_pos(k) = c*9;         ytick_lab{k} = sprintf('%s9', col_letter);
    end
    [ytick_pos, sort_idx] = sort(ytick_pos);
    ytick_lab = ytick_lab(sort_idx);

    fig = figure('Color', 'w', 'Position', [100 100 1100 700]);

    % --- Subplot 1: HD-sEMG (126 каналов) ---
    subplot(2,1,1);
    imagesc(time_vec, 1:n_hd, act_matrix);
    axis xy;
    colormap(jet);
    clim_hd = prctile(act_matrix(:), [1, 99]);
    if clim_hd(2) > clim_hd(1); caxis(clim_hd); end
    cb1 = colorbar; cb1.Label.String = colorbar_label;
    set(gca, 'YTick', ytick_pos, 'YTickLabel', ytick_lab, 'FontSize', 6);
    ylabel('HD-sEMG электрод');
    if nargin >= 4 && ~isempty(title_str)
        title(sprintf('HD-sEMG активность во времени (окно %d мс, перекрытие %.0f%%): %s', ...
            window_ms, overlap_ratio*100, title_str), 'Interpreter', 'none');
    else
        title(sprintf('HD-sEMG активность во времени (окно %d мс, перекрытие %.0f%%)', window_ms, overlap_ratio*100));
    end

    % --- Subplot 2: монополярные электроды EXG1-8 ---
    subplot(2,1,2);
    imagesc(time_vec, 1:n_exg, exg_matrix);
    axis xy;
    colormap(jet);
    clim_exg = prctile(exg_matrix(:), [1, 99]);
    if clim_exg(2) > clim_exg(1); caxis(clim_exg); end
    cb2 = colorbar; cb2.Label.String = colorbar_label;
    set(gca, 'YTick', 1:n_exg, 'YTickLabel', arrayfun(@(k) sprintf('EXG%d',k), 1:n_exg, 'UniformOutput', false));
    xlabel('Время, с'); ylabel('Монополярный электрод');
    title('Монополярные электроды (EXG1-EXG8)');

    fprintf('  [диагностика] HD: min=%.3f max=%.3f | EXG: min=%.3f max=%.3f | окон=%d (шаг=%d отсч.)\n', ...
        min(act_matrix(:)), max(act_matrix(:)), min(exg_matrix(:)), max(exg_matrix(:)), n_windows, step_samples);
end
