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

% --- Схема расположения электродов (справочная картинка, без данных) ---
fig_layout = plot_electrode_layout();
saveas(fig_layout, fullfile(figures_dir, 'fig_1_2_electrode_layout.png'));

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

% --- Расчёт RMS-сетки (без картинки) + heatmap во времени (окна 128мс) ---
% для 3 движений: кулак(04), pinch(08), index flexion(06)
movements_to_compare = struct( ...
    'code', {'04', '08', '06'}, ...
    'label', {'Кулак (fist)', 'Щипковый захват (pinch)', 'Сгибание пальца (index flexion)'});

heatmap_grids = cell(1, numel(movements_to_compare));
time_heatmaps = cell(1, numel(movements_to_compare));
exg_heatmaps = cell(1, numel(movements_to_compare));
d_all = cell(1, numel(movements_to_compare));   % сохраняем данные каждого движения для 1.3/1.4
window_ms = 128;   % длина окна для временного heatmap

for i = 1:numel(movements_to_compare)
    fname = sprintf('detop_exp01_subj01_Sess1_%s_01.mat', movements_to_compare(i).code);
    fpath = fullfile(data_dir, 'subj01', fname);

    d = load_data(fpath);
    d_all{i} = d;

    % Только расчёт RMS-сетки, без построения графика (нужно для 1.4/1.8)
    heatmap_grids{i} = calc_hdsemg_rms_grid(d.emg);

    % Heatmap во времени: HD-sEMG (126) + EXG (8) в одном окне, 2 subplot
    [act_matrix, exg_matrix, time_vec, fig_time] = plot_hdsemg_time_heatmap(d.emg, d.fs_emg, window_ms, movements_to_compare(i).label);
    time_heatmaps{i} = act_matrix;
    exg_heatmaps{i} = exg_matrix;
    saveas(fig_time, fullfile(figures_dir, sprintf('fig_1_2_time_heatmap_%s.png', movements_to_compare(i).code)));
end

fprintf('\nВременные heatmap-ы (окно %d мс) сохранены для движений: %s\n', ...
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
% По заданию: "для одного жеста" - берём index flexion (код 06), самый
% наглядный случай: один явный датчик пальца <-> один явный EMG-канал.

target_movement_idx = 3;   % movements_to_compare(3) = Сгибание пальца (index flexion)
di = d_all{target_movement_idx};
mov_label = movements_to_compare(target_movement_idx).label;

glove_names = cellstr(di.channels_glove);
t_glove_i = (0:size(di.glove, 2)-1) / di.fs_glove;
t_emg_i    = (0:size(di.emg, 2)-1)   / di.fs_emg;

% --- Все 18 каналов перчатки для этого жеста (одно окно) ---
figure('Color','w');
plot(t_glove_i, di.glove);
legend(glove_names, 'Location', 'eastoutside', 'FontSize', 7);
xlabel('Время, с'); ylabel('Нормализованный угол (0-1)');
title(sprintf('CyberGlove III: все 18 каналов, движение "%s"', mov_label));
saveas(gcf, fullfile(figures_dir, 'fig_1_4_glove_all_channels.png'));

% --- Топ EMG-канал (макс. RMS) для этого движения ---
grid_i = heatmap_grids{target_movement_idx};
[~, lin_idx] = max(grid_i(:));
[row_i, col_i] = ind2sub(size(grid_i), lin_idx);
emg_ch = (col_i - 1) * 9 + row_i;

% --- Подходящие датчики перчатки (Index_Inner, Index_Middle) ---
target_names = {'Index_Inner', 'Index_Middle'};
glove_idx = find(ismember(glove_names, target_names));

% --- Наложение: 1 EMG-канал + 1-2 датчика перчатки (одно окно) ---
figure('Color','w');
yyaxis left;
plot(t_emg_i, di.emg(emg_ch, :), 'b');
ylabel(sprintf('EMG: %s', strtrim(di.channels_emg(emg_ch, :))));

yyaxis right;
hold on;
colors = lines(numel(glove_idx));
for k = 1:numel(glove_idx)
    plot(t_glove_i, di.glove(glove_idx(k), :), 'Color', colors(k,:), 'LineWidth', 1.3);
end
ylabel('Перчатка: угол сустава (0-1)');

xlabel('Время, с');
legend([{sprintf('EMG %s (топ)', strtrim(di.channels_emg(emg_ch,:)))}, glove_names(glove_idx)'], ...
    'Location', 'best', 'FontSize', 8);
title(sprintf('ЭМГ (топ-канал) + перчатка: движение "%s"', mov_label));
saveas(gcf, fullfile(figures_dir, 'fig_1_4_emg_glove_overlay.png'));

fprintf('\n[1.4] %s:\n  EMG-канал (топ): %s\n  Датчики перчатки: %s\n', ...
    mov_label, strtrim(di.channels_emg(emg_ch, :)), strjoin(glove_names(glove_idx), ', '));

%% --- 1.5. Частотный анализ ---

fs = data.fs_emg;   % 2048 Гц (из п.1.1)

% --- PSD для 2 каналов (1 HD + 1 EXG), движение "Кулак" - подтвердить диапазон 10-500 Гц ---
di_freq = d_all{1};   % Кулак (fist)
ch_hd_freq  = 1;      % пример HD-канала (A1)
ch_exg_freq = 127;    % EXG1, первый монополярный

window_psd   = hamming(512);
noverlap_psd = 256;
nfft_psd     = 1024;

[pxx_hd, f_hd]   = pwelch(di_freq.emg(ch_hd_freq, :),  window_psd, noverlap_psd, nfft_psd, fs);
[pxx_exg, f_exg] = pwelch(di_freq.emg(ch_exg_freq, :), window_psd, noverlap_psd, nfft_psd, fs);

figure('Color','w');
plot(f_hd, 10*log10(pxx_hd), 'b'); hold on;
plot(f_exg, 10*log10(pxx_exg), 'r');
xline(10, '--k', '10 Гц');
xline(500, '--k', '500 Гц');
xlim([0, fs/2]);
xlabel('Частота, Гц'); ylabel('PSD, дБ/Гц');
legend({sprintf('HD-канал %s', strtrim(di_freq.channels_emg(ch_hd_freq,:))), ...
        sprintf('EXG-канал %s', strtrim(di_freq.channels_emg(ch_exg_freq,:)))});
title('PSD: диапазон значимой энергии сигнала (движение "Кулак")');
saveas(gcf, fullfile(figures_dir, 'fig_1_5_psd.png'));

fprintf('\n[1.5] PSD построена для каналов %s и %s\n', ...
    strtrim(di_freq.channels_emg(ch_hd_freq,:)), strtrim(di_freq.channels_emg(ch_exg_freq,:)));

% --- Медианная частота: сравнение по движениям (топ-канал каждого) ---
med_freqs = zeros(1, numel(movements_to_compare));
for i = 1:numel(movements_to_compare)
    di = d_all{i};
    grid_i = heatmap_grids{i};
    [~, lin_idx] = max(grid_i(:));
    [row_i, col_i] = ind2sub(size(grid_i), lin_idx);
    ch_i = (col_i - 1) * 9 + row_i;

    med_freqs(i) = medfreq(di.emg(ch_i, :), di.fs_emg);
end

figure('Color','w');
bar(med_freqs);
set(gca, 'XTickLabel', {movements_to_compare.label});
ylabel('Медианная частота, Гц');
title('Медианная частота (топ-канал) по движениям');
saveas(gcf, fullfile(figures_dir, 'fig_1_5_median_frequency.png'));

fprintf('\n[1.5] Медианные частоты (топ-канал каждого движения):\n');
for i = 1:numel(movements_to_compare)
    fprintf('  %s: %.1f Гц\n', movements_to_compare(i).label, med_freqs(i));
end
% --- АЧХ / амплитудный спектр одного канала ---
ch_amp = 1;
signal = double(di_freq.emg(ch_amp, :));
fs = di_freq.fs_emg;
signal = detrend(signal, 'constant');
N = length(signal);
Y = fft(signal);
P2 = abs(Y) / N;
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
f = fs * (0:floor(N/2)) / N;
figure('Color','w');
plot(f, P1, 'LineWidth', 1.2);
xlim([10 500]);
xlabel('Частота, Гц');
ylabel('Амплитуда');
title(sprintf('Амплитудный спектр ЭМГ — канал %s', ...
    strtrim(di_freq.channels_emg(ch_amp,:))));
grid on;
saveas(gcf, fullfile(figures_dir, 'fig_1_5_amplitude_spectrum.png'));

% --- Спектрограмма для одного канала (топ-канал "Кулак") ---
grid_fist = heatmap_grids{1};
[~, lin_idx] = max(grid_fist(:));
[row_f, col_f] = ind2sub(size(grid_fist), lin_idx);
ch_spectro = (col_f - 1) * 9 + row_f;

figure('Color','w');
spectrogram(di_freq.emg(ch_spectro, :), 256, 0, 2048, fs, 'yaxis');
title(sprintf('Спектрограмма: канал %s, движение "Кулак"', strtrim(di_freq.channels_emg(ch_spectro,:))));
saveas(gcf, fullfile(figures_dir, 'fig_1_5_spectrogram.png'));

%% --- Проверка артефактов в начале записи (15 случайных файлов) ---
rng(42);   % фиксируем seed для воспроизводимости выбора файлов
all_files = dir(fullfile(data_dir, 'subj01', 'detop_exp01_subj01_Sess*_*_*.mat'));
n_check = min(15, numel(all_files));
check_idx = randperm(numel(all_files), n_check);

check_window_ms = 50;   % проверяемое окно в самом начале записи, мс

fprintf('\n=== Проверка артефактов в начале записи (%d случайных файлов) ===\n', n_check);
fprintf('%-45s %10s %10s %8s\n', 'Файл', 'пик(нач)', 'медиана', 'отн-е');

figure('Color','w', 'Position', [50 50 1400 800]);
for k = 1:n_check
    fname = all_files(check_idx(k)).name;
    fpath = fullfile(all_files(check_idx(k)).folder, fname);
    dm = load_data(fpath);

    win_samples = round(check_window_ms / 1000 * dm.fs_emg);
    sig = mean(dm.emg(1:126, :).^2, 1);   % средняя мощность по всем HD-каналам

    first_win_max = max(sig(1:win_samples));
    rest_median   = median(sig(win_samples+1:end));
    ratio = first_win_max / (rest_median + eps);

    n_plot = min(numel(sig), round(dm.fs_emg * 1));   % первая секунда для графика
    t_plot = (0:n_plot-1) / dm.fs_emg;

    subplot(3, 5, k);
    plot(t_plot, sig(1:n_plot));
    xline(check_window_ms/1000, 'r--');
    title(sprintf('%s\nотн.=%.1fx', fname, ratio), 'Interpreter', 'none', 'FontSize', 6);
    xlabel('с', 'FontSize', 6); set(gca, 'FontSize', 5);

    flag = '';
    if ratio > 5
        flag = '  <-- ПОДОЗРИТЕЛЬНО';
    end
    fprintf('%-45s %10.4f %10.4f %7.2fx%s\n', fname, first_win_max, rest_median, ratio, flag);
end
sgtitle(sprintf('Проверка первых %d мс записи на артефакты (15 случ. файлов)', check_window_ms));
saveas(gcf, fullfile(figures_dir, 'fig_artifact_check.png'));

%% --- 1.6. Матрица корреляции ---
% Для трёх движений: кулак(04), pinch(08), index flexion(06)

% Физические координаты HD-каналов (row 1-9, col 1-14) - нужны для анализа
% "корреляция vs расстояние между электродами"
hd_row = zeros(1, 126); hd_col = zeros(1, 126);
for k = 1:126
    hd_col(k) = floor((k-1)/9) + 1;
    hd_row(k) = mod(k-1, 9) + 1;
end

corr_matrices = cell(1, numel(movements_to_compare));
suspect_channel_count = zeros(1, 126);   % счётчик "подозрительности" по каждому HD-каналу (по всем движениям)

for i = 1:numel(movements_to_compare)
    di = d_all{i};
    corr_matrix = corrcoef(di.emg');   % 134x134
    corr_matrices{i} = corr_matrix;

    figure('Color','w');
    imagesc(corr_matrix);
    axis image;
    colormap(jet);
    caxis([-1, 1]);
    cb = colorbar; cb.Label.String = 'Коэффициент корреляции';
    xlabel('Канал (1-126 HD, 127-134 EXG)');
    ylabel('Канал (1-126 HD, 127-134 EXG)');
    title(sprintf('Корреляция 134 каналов: движение "%s"', movements_to_compare(i).label));
    hold on;
    xline(126.5, 'k-', 'LineWidth', 1.5);
    yline(126.5, 'k-', 'LineWidth', 1.5);
    saveas(gcf, fullfile(figures_dir, sprintf('fig_1_6_correlation_%s.png', movements_to_compare(i).code)));
end

% --- Анализ: корреляция vs физическое расстояние между HD-электродами ---
% (отвечает на "соседние - ожидаемо, где перекрёстные помехи" из задания)
fprintf('\n=== 1.6: Корреляция в зависимости от расстояния между электродами ===\n');

crosstalk_thresh_dist = 5;    % "дальние" электроды: расстояние > 5 клеток сетки
crosstalk_thresh_corr = 0.7;  % подозрительно высокая корреляция для такой дистанции

for i = 1:numel(movements_to_compare)
    corr_hd = corr_matrices{i}(1:126, 1:126);

    pair_a = []; pair_b = []; distances = []; corrs = [];
    for a = 1:125
        for b = a+1:126
            d_ab = sqrt((hd_row(a)-hd_row(b))^2 + (hd_col(a)-hd_col(b))^2);
            pair_a(end+1) = a; pair_b(end+1) = b; %#ok<AGROW>
            distances(end+1) = d_ab; corrs(end+1) = corr_hd(a, b); %#ok<AGROW>
        end
    end

    neighbor_mask = distances <= 1.5;                       % соседние электроды
    mid_mask      = distances > 1.5 & distances <= 5;        % средняя дистанция
    far_mask      = distances > crosstalk_thresh_dist;       % дальние электроды

    fprintf('\n%s:\n', movements_to_compare(i).label);
    fprintf('  Соседние (d<=1.5):      средняя корр. = %.3f (n=%d)\n', mean(corrs(neighbor_mask)), sum(neighbor_mask));
    fprintf('  Средняя дистанция:      средняя корр. = %.3f (n=%d)\n', mean(corrs(mid_mask)), sum(mid_mask));
    fprintf('  Дальние (d>%.0f):        средняя корр. = %.3f (n=%d)\n', crosstalk_thresh_dist, mean(corrs(far_mask)), sum(far_mask));

    % --- Поиск подозрений на перекрёстные помехи: дальние пары с высокой корреляцией ---
    suspect_idx = find(far_mask & (corrs > crosstalk_thresh_corr));
    fprintf('  Подозрение на перекрёстные помехи (d>%.0f и corr>%.1f): %d пар\n', ...
        crosstalk_thresh_dist, crosstalk_thresh_corr, numel(suspect_idx));

    % Накопление счётчика "подозрительности" по каналам (для сводки после цикла)
    for s = 1:numel(suspect_idx)
        idx_pair = suspect_idx(s);
        suspect_channel_count(pair_a(idx_pair)) = suspect_channel_count(pair_a(idx_pair)) + 1;
        suspect_channel_count(pair_b(idx_pair)) = suspect_channel_count(pair_b(idx_pair)) + 1;
    end

    if ~isempty(suspect_idx)
        [~, sort_order] = sort(corrs(suspect_idx), 'descend');
        top_show = min(5, numel(suspect_idx));
        fprintf('  Топ-%d подозрительных пар:\n', top_show);
        for s = 1:top_show
            idx_pair = suspect_idx(sort_order(s));
            a = pair_a(idx_pair); b = pair_b(idx_pair);
            name_a = strtrim(di.channels_emg(a, :));
            name_b = strtrim(di.channels_emg(b, :));
            fprintf('    %s <-> %s: corr=%.3f, расстояние=%.1f\n', ...
                name_a, name_b, corrs(idx_pair), distances(idx_pair));
        end
    end
end

% --- Сводка: какие конкретные электроды чаще всего встречаются в подозрительных парах ---
% (суммарно по всем 3 движениям - если электрод неисправен, он будет "светиться" везде)
fprintf('\n=== 1.6: Сводный рейтинг электродов по частоте попадания в подозрительные пары ===\n');
[sorted_counts, sorted_ch] = sort(suspect_channel_count, 'descend');
top_suspect_n = 10;
for k = 1:top_suspect_n
    if sorted_counts(k) == 0
        break;
    end
    ch_idx = sorted_ch(k);
    col_letter = char('A' + hd_col(ch_idx) - 1);
    ch_name = sprintf('%s%d', col_letter, hd_row(ch_idx));
    fprintf('  %2d) канал %s (индекс %d): %d раз в подозрительных парах (по всем 3 движениям)\n', ...
        k, ch_name, ch_idx, sorted_counts(k));
end

%% --- 1.7. Нарезка каждого движения (сегментация) ---
% Все 13 движений, 1 повторение, subj01 - сетка 4x4

movement_codes_all  = arrayfun(@(x) sprintf('%02d', x), 1:13, 'UniformOutput', false);
movement_labels_all = {'Three-digit pinch','Cylinder grasp','Disc grasp','Fist', ...
    'Index+thumb trumpet','Index flexion','Middle+thumb trumpet','Pinch', ...
    'Thumb adduction','Thumb extension','Thumb flexion','Point','MRL flexion'};

segmentation_results = struct('code', cell(1,13), 'label', cell(1,13), ...
    'onset_s', cell(1,13), 'offset_s', cell(1,13), 'action_duration_s', cell(1,13));

plot_channels_1_7 = [1, 63];   % A1 (край сетки) + канал ближе к центру, для наглядности

figure('Color','w', 'Position', [50 50 1400 900]);
for m = 1:13
    fname = sprintf('detop_exp01_subj01_Sess1_%s_01.mat', movement_codes_all{m});
    fpath = fullfile(data_dir, 'subj01', fname);

    dm = load_data(fpath);
    [onset_t, offset_t] = detect_movement_onset(dm.emg, dm.fs_emg);
    t_dm = (0:size(dm.emg, 2)-1) / dm.fs_emg;

    subplot(4, 4, m);
    hold on;
    for ch = plot_channels_1_7
        plot(t_dm, dm.emg(ch, :));
    end
    if ~isnan(onset_t)
        xline(onset_t, 'g--', 'LineWidth', 1.2);
        xline(offset_t, 'r--', 'LineWidth', 1.2);
    end
    title(sprintf('%s (%s)', movement_labels_all{m}, movement_codes_all{m}), 'FontSize', 7);
    xlabel('с', 'FontSize', 6); set(gca, 'FontSize', 5);

    segmentation_results(m).code = movement_codes_all{m};
    segmentation_results(m).label = movement_labels_all{m};
    segmentation_results(m).onset_s = onset_t;
    segmentation_results(m).offset_s = offset_t;
    segmentation_results(m).action_duration_s = offset_t - onset_t;
end
sgtitle('Сегментация: покой -> действие -> покой (subj01, rep 1, 13 движений)');
saveas(gcf, fullfile(figures_dir, 'fig_1_7_segmentation_grid.png'));

fprintf('\n=== 1.7: Длительности движений (повторение 1) ===\n');
fprintf('%-30s %8s %8s %10s\n', 'Движение', 'начало', 'конец', 'длит-ть');
for m = 1:13
    fprintf('%-30s %7.2fс %7.2fс %8.2fс\n', segmentation_results(m).label, ...
        segmentation_results(m).onset_s, segmentation_results(m).offset_s, ...
        segmentation_results(m).action_duration_s);
end

durations_rep1 = [segmentation_results.action_duration_s];
fprintf('\n[1.7] Среднее=%.2fс, std=%.2fс, min=%.2fс, max=%.2fс\n', ...
    mean(durations_rep1, 'omitnan'), std(durations_rep1, 'omitnan'), ...
    min(durations_rep1), max(durations_rep1));

%% --- 1.7 (доп). Сегментация ВСЕХ файлов (13 движений x 6 повторений) ---
% Сводная "матрица": по вертикали - файл, по горизонтали - время,
% серый = вся запись (покой+действие+покой), синий = только действие,
% зелёная/красная точка = начало/конец действия

n_reps = 6;
all_segments = struct('code', {}, 'rep', {}, 'label', {}, 'onset_s', {}, 'offset_s', {}, 'total_duration_s', {});

for m = 1:13
    for r = 1:n_reps
        fname = sprintf('detop_exp01_subj01_Sess1_%s_%02d.mat', movement_codes_all{m}, r);
        fpath = fullfile(data_dir, 'subj01', fname);
        if ~isfile(fpath)
            continue;
        end
        dm = load_data(fpath);
        [onset_t, offset_t] = detect_movement_onset(dm.emg, dm.fs_emg);
        total_dur = size(dm.emg, 2) / dm.fs_emg;

        idx = numel(all_segments) + 1;
        all_segments(idx).code = movement_codes_all{m};
        all_segments(idx).rep = r;
        all_segments(idx).label = movement_labels_all{m};
        all_segments(idx).onset_s = onset_t;
        all_segments(idx).offset_s = offset_t;
        all_segments(idx).total_duration_s = total_dur;
    end
end

n_all = numel(all_segments);
figure('Color','w', 'Position', [50 50 1000 1400]);
hold on;
y_labels = cell(1, n_all);
for k = 1:n_all
    s = all_segments(k);
    y_labels{k} = sprintf('%s #%d', s.label, s.rep);

    plot([0, s.total_duration_s], [k, k], 'Color', [0.8 0.8 0.8], 'LineWidth', 1);   % вся запись
    if ~isnan(s.onset_s)
        plot([s.onset_s, s.offset_s], [k, k], 'b-', 'LineWidth', 4);   % активная фаза
        plot(s.onset_s, k, 'g.', 'MarkerSize', 10);
        plot(s.offset_s, k, 'r.', 'MarkerSize', 10);
    end
end
set(gca, 'YTick', 1:n_all, 'YTickLabel', y_labels, 'FontSize', 5, 'YDir', 'reverse');
xlabel('Время, с');
title('Сегментация всех файлов (13 движений x 6 повторений)');
ylim([0, n_all+1]);
saveas(gcf, fullfile(figures_dir, 'fig_1_7_all_segments_matrix.png'));

fprintf('\n[1.7 доп] Обработано файлов: %d\n', n_all);

fprintf('\n=== 1.7 (доп): Разброс длительности действия по повторениям ===\n');
for m = 1:13
    mask = strcmp({all_segments.code}, movement_codes_all{m});
    durs = [all_segments(mask).offset_s] - [all_segments(mask).onset_s];
    fprintf('%-30s: среднее=%.2fс, std=%.2fс, min=%.2fс, max=%.2fс (n=%d)\n', ...
        movement_labels_all{m}, mean(durs, 'omitnan'), std(durs, 'omitnan'), ...
        min(durs), max(durs), sum(mask));
end

%% --- 1.8. Сравнительная таблица по 13 движениям ---
n_reps_18 = 6;
comparison_data = struct('code', cell(1,13), 'label', cell(1,13), ...
    'rms_uV', cell(1,13), 'mean_amp_uV', cell(1,13), 'median_freq_Hz', cell(1,13), 'duration_s', cell(1,13));

for m = 1:13
    rms_vals      = nan(1, n_reps_18);
    mean_amp_vals = nan(1, n_reps_18);
    medfreq_vals  = nan(1, n_reps_18);
    dur_vals      = nan(1, n_reps_18);

    for r = 1:n_reps_18
        fname = sprintf('detop_exp01_subj01_Sess1_%s_%02d.mat', movement_codes_all{m}, r);
        fpath = fullfile(data_dir, 'subj01', fname);
        if ~isfile(fpath)
            continue;
        end
        dm = load_data(fpath);

        rms_vals(r)      = sqrt(mean(dm.emg(:).^2));
        mean_amp_vals(r) = mean(abs(dm.emg(:)));

        avg_signal = mean(dm.emg, 1);
        medfreq_vals(r) = medfreq(avg_signal, dm.fs_emg);

        [onset_t, offset_t] = detect_movement_onset(dm.emg, dm.fs_emg);
        dur_vals(r) = offset_t - onset_t;
    end

    comparison_data(m).code           = movement_codes_all{m};
    comparison_data(m).label          = movement_labels_all{m};
    comparison_data(m).rms_uV         = mean(rms_vals, 'omitnan');
    comparison_data(m).mean_amp_uV    = mean(mean_amp_vals, 'omitnan');
    comparison_data(m).median_freq_Hz = mean(medfreq_vals, 'omitnan');
    comparison_data(m).duration_s     = mean(dur_vals, 'omitnan');
end

table_1_8 = struct2table(comparison_data);
disp(table_1_8);
writetable(table_1_8, fullfile(results_dir, 'table_1_8_movements_comparison.csv'));

% --- Гистограммы по каждому столбцу ---
figure('Color','w');
subplot(2,2,1); histogram(table_1_8.rms_uV, 6);          title('RMS амплитуда'); xlabel('мкВ');
subplot(2,2,2); histogram(table_1_8.mean_amp_uV, 6);     title('Средняя амплитуда'); xlabel('мкВ');
subplot(2,2,3); histogram(table_1_8.median_freq_Hz, 6);  title('Медианная частота'); xlabel('Гц');
subplot(2,2,4); histogram(table_1_8.duration_s, 6);      title('Длительность движения'); xlabel('с');
sgtitle('Распределение характеристик по 13 движениям');
saveas(gcf, fullfile(figures_dir, 'fig_1_8_histograms.png'));

fprintf('\n[1.8] Таблица сохранена: results/table_1_8_movements_comparison.csv\n');

%% --- 1.9. Итоговый вывод ---
% TODO
