function [act_matrix, time_vec, fig] = plot_hdsemg_time_heatmap(emg, fs_emg, window_ms, title_str)
% PLOT_HDSEMG_TIME_HEATMAP  Heatmap "канал x время" по HD-sEMG матрице:
%   для каждого канала считает RMS-амплитуду в неперекрывающихся окнах
%   длиной window_ms (мс), показывает изменение активности во времени.
%
%   [act_matrix, time_vec, fig] = PLOT_HDSEMG_TIME_HEATMAP(emg, fs_emg, window_ms, title_str)
%
%   Вход:
%     emg        - 134 x N (или 126 x N) массив EMG-сигналов
%     fs_emg     - частота дискретизации EMG, Гц (для SEEDS = 2048)
%     window_ms  - длина окна в мс (по заданию = 128)
%     title_str  - заголовок графика (например, название движения)
%
%   Выход:
%     act_matrix - 126 x n_windows матрица RMS-амплитуды (канал x окно)
%     time_vec   - 1 x n_windows вектор времени (центр каждого окна), с
%     fig        - handle фигуры

    if nargin < 3 || isempty(window_ms)
        window_ms = 128;
    end

    n_hd = 126;
    emg_hd = emg(1:n_hd, :);

    win_samples = round(window_ms / 1000 * fs_emg);   % 128 мс * 2048 Гц = 262 отсчёта
    n_samples   = size(emg_hd, 2);
    n_windows   = floor(n_samples / win_samples);

    act_matrix = zeros(n_hd, n_windows);
    for w = 1:n_windows
        idx = (w-1)*win_samples + 1 : w*win_samples;
        act_matrix(:, w) = sqrt(mean(emg_hd(:, idx).^2, 2));
    end

    % Время в секундах: центр каждого окна
    time_vec = ((1:n_windows) - 0.5) * win_samples / fs_emg;

    % --- Heatmap ---
    fig = figure('Color', 'w');
    imagesc(time_vec, 1:n_hd, act_matrix);
    axis xy;   % канал 1 внизу, растёт вверх (как в обычных EMG-heatmap)
    colormap(jet);

    % Робастные пределы цвета: обрезаем по перцентилям, чтобы редкие
    % выбросы (например, краевой артефакт фильтрации в первом окне)
    % не "съедали" динамический диапазон реальной мышечной активности.
    clim_vals = prctile(act_matrix(:), [1, 99]);
    if clim_vals(2) > clim_vals(1)
        caxis(clim_vals);
    end

    cb = colorbar;
    cb.Label.String = 'RMS амплитуда';

    xlabel('Время, с');
    ylabel('Канал (индекс 1-126)');

    % Диагностика выбросов (особенно первое/последнее окно —
    % типичное место краевых артефактов офлайн-фильтрации)
    fprintf('  [диагностика] act_matrix: min=%.3f, max=%.3f, окно#1 max=%.3f, окно#%d(посл.) max=%.3f\n', ...
        min(act_matrix(:)), max(act_matrix(:)), ...
        max(act_matrix(:,1)), n_windows, max(act_matrix(:,end)));

    if nargin >= 4 && ~isempty(title_str)
        title(sprintf('HD-sEMG активность во времени (окно %d мс): %s', window_ms, title_str), ...
            'Interpreter', 'none');
    else
        title(sprintf('HD-sEMG активность во времени (окно %d мс)', window_ms));
    end
end
