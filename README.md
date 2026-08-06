# HD-sEMG Analysis — EDA набора SEEDS

Разведочный анализ (EDA) датасета SEEDS: синхронные записи HD-sEMG (126 каналов, сетка 9×14),
8 монополярных электродов-разгибателей и кинематики (CyberGlove III, 18 каналов) при выполнении
13 жестов кисти.

Источник датасета: Matran-Fernandez et al., *"SEEDS: a large-scale multimodal dataset for
human-computer interaction"*, Scientific Data (2019). https://www.nature.com/articles/s41597-019-0200-9

## Структура проекта

```
hdsemg-analysis/
├── data/
│   └── subj01/             # .mat файлы испытуемого (не коммитятся, см. .gitignore)
│       README.md           # ссылка на источник данных
├── code/
│   ├── eda_seeds.m         # основной EDA-скрипт (пункты 1.1-1.9)
│   └── utils/
│       ├── load_data.m                # загрузка .mat, вывод структуры переменных (1.1)
│       ├── plot_hdsemg_heatmap.m      # heatmap 9x14, RMS за всё движение (1.2)
│       └── plot_hdsemg_time_heatmap.m # heatmap канал x время, окно 128мс (1.2)
├── docs/
│   ├── report.docx         # итоговый отчёт
│   └── figures/            # сохранённые графики (генерируются скриптом)
├── results/                # таблицы для отчёта (генерируются скриптом)
├── .gitignore
└── README.md
```

## Как запустить

1. Скачать данные SEEDS (см. `data/README.md`), распаковать в `data/subj01/`
2. Открыть MATLAB, в редакторе открыть `code/eda_seeds.m`
3. Запустить (F5) — скрипт сам добавит `code/utils` в путь и создаст `results/`, `docs/figures/`,
   если их нет
4. Графики появятся в `docs/figures/`, таблицы — в `results/`

## Данные

Файл: `detop_exp01_subj{NN}_Sess{S}_{MM}_{RR}.mat`, где NN — испытуемый, S — сессия (1-3),
MM — движение (01-13), RR — повторение (1-6).

Переменные внутри: `emg` (134×N: 126 HD-sEMG + 8 EXG), `glove` (18×N), `fs_emg` (2048 Гц),
`fs_glove` (256 Гц), `channels_emg`, `channels_glove`, `subject`, `session`, `date`, `movement`, `speed`.
