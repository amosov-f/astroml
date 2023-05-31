import io
import json
import logging
import os
import re
from gzip import GzipFile
from urllib.request import urlopen
import yt.wrapper as yt

import pandas as pd
from numpy import float64, int64
from yt.wrapper import Format, DsvFormat, SchemafulDsvFormat

logging.basicConfig(level=logging.DEBUG)

GAIA_URL = 'http://cdn.gea.esac.esa.int/Gaia/gdr2/vari_cepheid/csv/' #'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/'

yt.config['proxy']['url'] = 'hahn'

def get_gaia_urls():
    # return [
    #     'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/GaiaSource_141120118403196800_141350912766009216.csv.gz',
    #     'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/GaiaSource_6105976231005817088_6106080895067888256.csv.gz'
    # ]
    urls = list(map(lambda filename: GAIA_URL + filename, get_gaia_names()))
    logging.info('Total ' + str(len(urls)) + ' urls')
    return urls


def get_gaia_names():
    string = urlopen(GAIA_URL).read().decode('utf-8')
    pattern = re.compile('VariCepheid.*\.csv\.gz"')# re.compile('GaiaSource.*\.csv\.gz"')  # the pattern actually creates duplicates in the list
    filelist = pattern.findall(string)
    return list(map(lambda filename: filename[:-1], filelist))

columns = [
    'source_id',
    'ra',
    'ra_error',
    'dec',
    'dec_error',
    'parallax',
    'parallax_error',
    'pmra',
    'pmra_error',
    'pmdec',
    'pmdec_error',
    'phot_g_mean_mag',
    'phot_rp_mean_mag',
    'radial_velocity',
    'radial_velocity_error'
]

def process(url, i):
    logging.info('Processing ' + str(i) + ': ' + url)
    table = '//home/cocainum/amosov-f/gaia/vari_cepheid'
    df = pd.read_csv(url, compression='gzip')
    if i == 0:
        schema = []
        for name, type in zip(df.dtypes.keys(), df.dtypes):
            # if name in columns:
            schema.append({
                'name': name,
                'type': to_yt_type(name, type)
            })
        yt.create('table', table, recursive=True, force=True, attributes={
            'schema': schema
        })
    file = 'filtered.jsonl'
    # df = df[columns]
    df.to_json(file, orient='records', lines=True)
    yt.write_table(yt.TablePath(table, append=True), open(file, 'rb'), raw=True, format='json')


def to_yt_type(name, type):
    if name == 'multi_mode_best_classification':
        return 'string'
    if type == float64:
        return 'double'
    if type == int64:
        return 'int64'
    if type == bool:
        return 'boolean'
    return 'string'

def main():
    problem_urls_file = 'problem_urls.txt'
    f = open(problem_urls_file, 'w')
    urls = get_gaia_urls()
    for i in range(0, len(urls)):
        url = urls[i]
        try:
            process(url, i)
        except:
            logging.exception('Can\'t process dataset #' + str(i) + ': ' + url)
            f.write(str(i) + ' ' + url + '\n')
            f.flush()
        print()


    # print(get_gaia_urls())




if __name__ == '__main__':
    main()
