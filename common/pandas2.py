import time

import pandas


def split(dataset):
    X = dataset.loc[:, dataset.columns != 'y']
    y = dataset.y
    return X, y


def series(**kwargs):
    return kwargs


def fast_apply(df: pandas.DataFrame, func):
    start = time.time()
    new_columns = {}
    for _, row in df.iterrows():
        new_values = func(row)
        for k, v in new_values.items():
            if k not in new_columns:
                new_columns[k] = []
            new_columns[k].append(v)
    print(f'Fast apply {func} finished in {int(time.time() - start)}s')
    return pandas.DataFrame(new_columns)
