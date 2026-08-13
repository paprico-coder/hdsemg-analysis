function [grid_rms, fig] = plot_hdsemg_heatmap(emg, channels_emg, title_str)
%   PLOT_HDSEMG_HEATMAP  Строит heatmap 9x14 RMS-амплитуды HD-sEMG матрицы.
%
%   [grid_rms, fig] = PLOT_HDSEMG_HEATMAP(emg, channels_emg, title_str)
%
%   Вход:
%     emg          - 134 x N (или 126 x N) массив EMG-сигналов
%     channels_emg - имена каналов (134x4 или 126x4 char array), нужны
%                    только чтобы вычислить (row, col) по имени вида 'MA1'
%     title_str    - заголовок графика (например, название движения)
%
%   Выход:
%     grid_rms - 9x14 матрица RMS-амплитуды (строка=1..9, столбец=A..N)
%     fig      - handle фигуры

    % Берём только первые 126 каналов (HD-sEMG матрица, без EXG1-8)
    n_hd = 126;
    emg_hd = emg(1:n_hd, :);

    % --- RMS по каждому из 126 каналов ---
    rms_vals = sqrt(mean(emg_hd.^2, 2));   % 126 x 1

    % --- Раскладываем 126 значений в сетку 9x14 ---
    % Порядок каналов: MA1..MA9 (col=A), MB1..MB9 (col=B), ..., MN1..MN9 (col=N)
    % col = floor((k-1)/9)+1 (1..14 -> A..N), row = mod(k-1,9)+1 (1..9)
    grid_rms = nan(9, 14);
    for k = 1:n_hd
        col = floor((k-1)/9) + 1;
        row = mod(k-1, 9) + 1;
        grid_rms(row, col) = rms_vals(k);
    end

    % --- Heatmap ---
    fig = figure('Color', 'w');
    imagesc(grid_rms);
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'RMS амплитуда';
    axis image;

    col_labels = arrayfun(@(c) char('A' + c - 1), 1:14, 'UniformOutput', false);
    set(gca, 'XTick', 1:14, 'XTickLabel', col_labels, ...
             'YTick', 1:9,  'YTickLabel', 1:9);
    xlabel('Столбец (A-N)');
    ylabel('Строка (1-9)');

    if nargin >= 3 && ~isempty(title_str)
        title(sprintf('HD-sEMG RMS heatmap: %s', title_str), 'Interpreter', 'none');
    else
        title('HD-sEMG RMS heatmap');
    end
end
