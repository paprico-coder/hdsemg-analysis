%% eda_seeds.m
% Основной скрипт EDA (разведочного анализа) для датасета SEEDS.
% Пошагово покрывает пункты 1.1-1.9 задания.
%
% Запускать из корня репозитория hdsemg-analysis/ (или поправить addpath ниже).

clear; clc; close all;

%% --- Настройка путей ---
% Определяем корень репозитория относительно расположения этого файла,
% а не текущей папки MATLAB (Current Folder может быть где угодно).
script_dir = fileparts(mfilename('fullpath'));   % .../hdsemg-analysis/code
repo_root  = fileparts(script_dir);              % .../hdsemg-analysis

addpath(fullfile(script_dir, 'utils'));

data_dir    = fullfile(repo_root, 'data');
results_dir = fullfile(repo_root, 'results');
figures_dir = fullfile(repo_root, 'docs', 'figures');

if ~isfolder(results_dir); mkdir(results_dir); end
if ~isfolder(figures_dir); mkdir(figures_dir); end

%% --- 1.1. Загрузка и структура .mat файла ---
% Берём один файл для разбора структуры (subject 01, session 1, movement 01, rep 01)
% Файлы лежат в подпапке по номеру испытуемого: data/subj01/...
example_file = fullfile(data_dir, 'subj01', 'detop_exp01_subj01_Sess1_01_01.mat');

data = load_data(example_file);

% Таблица переменных (Variable/Size/Type/Description) уже выведена в
% Command Window внутри load_data.m и сохранена в base workspace как
% seeds_var_table. Сохраняем её в results/ для отчёта:
writetable(seeds_var_table, fullfile(results_dir, 'table_1_1_variables.csv'));

%% --- 1.2. HD-sEMG матрица (126 каналов, 9x14) ---

% --- Временные ряды 6-9 каналов из разных зон матрицы (на файле из 1.1) ---
t_example = (0:size(data.emg, 2)-1) / data.fs_emg;   % ось времени, с

figure('Color','w');
sample_channels = [1, 9, 63, 71, 118, 126];  % углы сетки (A1,A9) + центр + угол (N1,N9)
hold on;
for ch = sample_channels
    plot(t_example, data.emg(ch, :));
end
legend(cellstr(data.channels_emg(sample_channels, :)));
xlabel('Время, с'); ylabel('Амплитуда');
title('HD-sEMG: примеры каналов из разных зон матрицы');
saveas(gcf, fullfile(figures_dir, 'fig_1_2_channels_timeseries.png'));

% --- Heatmap (RMS за всё движение) + heatmap во времени (окна 128мс) ---
% для 3 движений: кулак(04), pinch(08), index flexion(06)
movements_to_compare = struct( ...
    'code', {'04', '08', '06'}, ...
    'label', {'Кулак (fist)', 'Щипковый захват (pinch)', 'Сгибание пальца (index flexion)'});

heatmap_grids = cell(1, numel(movements_to_compare));
time_heatmaps = cell(1, numel(movements_to_compare));
window_ms = 128;   % длина окна для временного heatmap

for i = 1:numel(movements_to_compare)
    fname = sprintf('detop_exp01_subj01_Sess1_%s_01.mat', movements_to_compare(i).code);
    fpath = fullfile(data_dir, 'subj01', fname);

    d = load_data(fpath);

    % Статический heatmap (одно значение RMS на канал за всю запись)
    [grid_rms, fig_hm] = plot_hdsemg_heatmap(d.emg, d.channels_emg, movements_to_compare(i).label);
    heatmap_grids{i} = grid_rms;
    saveas(fig_hm, fullfile(figures_dir, sprintf('fig_1_2_heatmap_%s.png', movements_to_compare(i).code)));

    % Heatmap во времени (канал x время, окно 128 мс)
    [act_matrix, time_vec, fig_time] = plot_hdsemg_time_heatmap(d.emg, d.fs_emg, window_ms, movements_to_compare(i).label);
    time_heatmaps{i} = act_matrix;
    saveas(fig_time, fullfile(figures_dir, sprintf('fig_1_2_time_heatmap_%s.png', movements_to_compare(i).code)));
end

fprintf('\nHeatmap-ы (статич. + временные, окно %d мс) сохранены для движений: %s\n', ...
    window_ms, strjoin({movements_to_compare.label}, ', '));

% --- Анализ для вывода по 1.2: топ-каналы и сравнение зон активации ---
fprintf('\n=== Анализ активных зон (для вывода п.1.2) ===\n');

top_n = 5;
for i = 1:numel(movements_to_compare)
    grid = heatmap_grids{i};                 % 9x14
    [sorted_vals, lin_idx] = sort(grid(:), 'descend');
    [rows, cols] = ind2sub(size(grid), lin_idx(1:top_n));

    fprintf('\n%s:\n', movements_to_compare(i).label);
    for k = 1:top_n
        col_letter = char('A' + cols(k) - 1);
        fprintf('  %d) канал %s%d (строка %d, столбец %s), RMS = %.4f\n', ...
            k, col_letter, rows(k), rows(k), col_letter, sorted_vals(k));
    end
end

% --- Попарная корреляция сеток активации между движениями ---
fprintf('\n=== Схожесть зон активации между движениями (корреляция сеток) ===\n');
n_mov = numel(movements_to_compare);
for i = 1:n_mov
    for j = i+1:n_mov
        r = corr(heatmap_grids{i}(:), heatmap_grids{j}(:));
        fprintf('  %s  vs  %s:  r = %.3f\n', ...
            movements_to_compare(i).label, movements_to_compare(j).label, r);
    end
end

%% --- 1.3. Монополярные электроды (8 каналов) ---

% --- Все 8 EXG-каналов на одном графике (файл из 1.1) ---
figure('Color','w');
exg_idx = 127:134;   % последние 8 каналов = EXG1..EXG8
hold on;
for ch = exg_idx
    plot(t_example, data.emg(ch, :));
end
legend(cellstr(data.channels_emg(exg_idx, :)));
xlabel('Время, с'); ylabel('Амплитуда');
title('Монополярные электроды (EXG1-EXG8, разгибатели)');
saveas(gcf, fullfile(figures_dir, 'fig_1_3_exg_channels.png'));

% --- Сравнение с HD-матрицей: сгибатель (топ-канал из HD) vs разгибатель (EXG) ---
% Берём топ-1 HD-канал по RMS (из анализа heatmap_grids движения "Кулак")
grid_fist = heatmap_grids{1};   % соответствует movements_to_compare(1) = Кулак
[~, lin_idx] = max(grid_fist(:));
[top_row, top_col] = ind2sub(size(grid_fist), lin_idx);
top_hd_channel = (top_col - 1) * 9 + top_row;   % индекс канала 1..126

% RMS по EXG-каналам для того же движения (кулак)
d_fist = load_data(fullfile(data_dir, 'subj01', 'detop_exp01_subj01_Sess1_04_01.mat'));
exg_rms = sqrt(mean(d_fist.emg(exg_idx, :).^2, 2));
[~, top_exg_local] = max(exg_rms);
top_exg_channel = exg_idx(top_exg_local);

t_fist = (0:size(d_fist.emg, 2)-1) / d_fist.fs_emg;

figure('Color','w');
yyaxis left;
plot(t_fist, d_fist.emg(top_hd_channel, :), 'b');
ylabel(sprintf('Амплитуда: %s (сгибатель)', strtrim(d_fist.channels_emg(top_hd_channel, :))));

yyaxis right;
plot(t_fist, d_fist.emg(top_exg_channel, :), 'r');
ylabel(sprintf('Амплитуда: %s (разгибатель)', strtrim(d_fist.channels_emg(top_exg_channel, :))));

xlabel('Время, с');
legend({'HD-канал (сгибатель)', 'EXG-канал (разгибатель)'}, 'Location', 'best');
title('Сравнение сгибатель (HD) vs разгибатель (EXG): движение "Кулак"');
saveas(gcf, fullfile(figures_dir, 'fig_1_3_flexor_vs_extensor.png'));

fprintf('\n[1.3] Топ HD-канал (сгибатель) для "Кулак": %s, топ EXG-канал (разгибатель): %s\n', ...
    strtrim(d_fist.channels_emg(top_hd_channel, :)), strtrim(d_fist.channels_emg(top_exg_channel, :)));

%% --- 1.4. CyberGlove III (18 каналов) ---
% Используем последний загруженный в цикле 1.2/1.3 файл: движение "Сгибание пальца" (код 06)
% В переменной d уже лежат данные этого движения.

t_glove = (0:size(d.glove, 2)-1) / d.fs_glove;   % ось времени, с

% --- Все 18 каналов перчатки на одном графике ---
figure('Color','w');
plot(t_glove, d.glove);
glove_names = cellstr(d.channels_glove);
legend(glove_names, 'Location', 'eastoutside', 'FontSize', 7);
xlabel('Время, с'); ylabel('Нормализованный угол (0-1)');
title(sprintf('CyberGlove III: все 18 каналов, движение "%s"', movements_to_compare(3).label));
saveas(gcf, fullfile(figures_dir, 'fig_1_4_glove_all_channels.png'));

fprintf('\n[1.4] Каналы перчатки (channels_glove):\n');
disp(glove_names);

% --- Автопоиск канала перчатки, соответствующего указательному пальцу ---
target_keyword = 'Index';
idx_glove_matches = find(contains(glove_names, target_keyword, 'IgnoreCase', true));

if isempty(idx_glove_matches)
    warning('Не найден канал перчатки по ключу "%s" - проверьте channels_glove вручную и укажите индекс.', target_keyword);
else
    fprintf('\n Найдены каналы перчатки по ключу "%s":\n', target_keyword);
    for k = 1:numel(idx_glove_matches)
        fprintf('  idx=%d: %s\n', idx_glove_matches(k), glove_names{idx_glove_matches(k)});
    end

    glove_ch = idx_glove_matches(1);   % берём первый найденный (обычно MCP-сустав)

    % Соответствующий EMG-канал: топ HD-канал для движения "Сгибание пальца" (heatmap_grids{3})
    grid_idxflex = heatmap_grids{3};
    [~, lin_idx] = max(grid_idxflex(:));
    [row_i, col_i] = ind2sub(size(grid_idxflex), lin_idx);
    emg_ch = (col_i - 1) * 9 + row_i;

    % --- Наложение: EMG-канал + датчик сгибания пальца перчатки ---
    t_d = (0:size(d.emg, 2)-1) / d.fs_emg;   % ось времени для EMG этого файла

    figure('Color','w');
    yyaxis left;
    plot(t_d, d.emg(emg_ch, :), 'b');
    ylabel(sprintf('EMG: %s', strtrim(d.channels_emg(emg_ch, :))));

    yyaxis right;
    plot(t_glove, d.glove(glove_ch, :), 'r', 'LineWidth', 1.5);
    ylabel(sprintf('Перчатка: %s', glove_names{glove_ch}));

    xlabel('Время, с');
    legend({'EMG (сгибатель)', 'Перчатка (угол сустава)'}, 'Location', 'best');
    title(sprintf('Наложение ЭМГ + перчатка: движение "%s"', movements_to_compare(3).label));
    saveas(gcf, fullfile(figures_dir, 'fig_1_4_emg_glove_overlay.png'));

    fprintf('\n[1.4] Наложены: EMG-канал %s vs перчатка %s\n', ...
        strtrim(d.channels_emg(emg_ch, :)), glove_names{glove_ch});
end

%% --- 1.5. Частотный анализ ---
% TODO

%% --- 1.6. Матрица корреляции ---
% TODO

%% --- 1.7. Нарезка каждого движения (сегментация) ---
% TODO

%% --- 1.8. Сравнительная таблица по 13 движениям ---
% TODO

%% --- 1.9. Итоговый вывод ---
% TODO
