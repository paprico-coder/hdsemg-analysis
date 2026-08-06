function data = load_data(filepath)
%   LOAD_DATA  Загружает один .mat файл датасета SEEDS и выводит
%   структуру всех переменных (размерность, тип, описание) в виде таблицы.
%
%   data = LOAD_DATA(filepath) — загружает файл по пути filepath
%   Пример:
%       data = load_data('data/detop_exp01_subj01_Sess1_01_01.mat');
%       data.emg     % 134 x N матрица EMG
%       data.glove   % 18 x N матрица перчатки

    if ~isfile(filepath)
        error('load_data:fileNotFound', 'Файл не найден: %s', filepath);
    end

    % Загружаем всё содержимое .mat в struct (не в workspace!)
    data = load(filepath);

    % --- Словарь описаний по документации SEEDS (Matran-Fernandez et al., 2019) ---
    descriptions = struct( ...
        'subject',       'номер испытуемого', ...
        'session',       'номер сессии (1-3)', ...
        'date',          'дата записи, формат YYYYMMDD', ...
        'movement',      'название выполненного движения', ...
        'fs_emg',        'частота дискретизации EMG, Гц (2048)', ...
        'fs_glove',      'частота дискретизации перчатки, Гц (256)', ...
        'speed',         'скорость выполнения (fast/slow)', ...
        'emg',           'EMG-сигналы: 134 канала x N отсчётов (126 HD-sEMG + 8 моно)', ...
        'channels_emg',  'имена каналов EMG по порядку строк emg', ...
        'glove',         'сигналы перчатки: 18 каналов x N отсчётов, нормализовано 0-1', ...
        'channels_glove','имена каналов перчатки по порядку строк glove' ...
    );

    % --- Собираем таблицу: имя, размерность, тип, описание ---
    fields = fieldnames(data);
    nVars = numel(fields);

    varNames   = cell(nVars, 1);
    varSizes   = cell(nVars, 1);
    varTypes   = cell(nVars, 1);
    varDescr   = cell(nVars, 1);

    for i = 1:nVars
        name = fields{i};
        val = data.(name);

        varNames{i} = name;
        varSizes{i} = mat2str(size(val));
        varTypes{i} = class(val);

        if isfield(descriptions, name)
            varDescr{i} = descriptions.(name);
        else
            varDescr{i} = '(нет в документации SEEDS)';
        end
    end

    T = table(varNames, varSizes, varTypes, varDescr, ...
        'VariableNames', {'Variable', 'Size', 'Type', 'Description'});

    fprintf('\n=== Структура файла: %s ===\n', filepath);
    disp(T);

    % Сохраняем таблицу для использования в отчёте (docx/results)
    assignin('base', 'seeds_var_table', T);
end
