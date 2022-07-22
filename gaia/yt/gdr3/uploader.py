import logging
import re
from urllib.request import urlopen

import pandas as pd

logging.basicConfig(level=logging.DEBUG)

GAIA_URL = 'http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/'


def get_gaia_urls():
    # return [
    #     'file:///Users/amosov-f/Documents/univer/astroml/gaia/yt/gdr3/GaiaSource_000000-003111.csv.gz'
    # ]
    # return [
    #     'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/GaiaSource_141120118403196800_141350912766009216.csv.gz',
    #     'http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/GaiaSource_6105976231005817088_6106080895067888256.csv.gz'
    # ]
    urls = list(map(lambda filename: GAIA_URL + filename, get_gaia_names()))
    logging.info('Total ' + str(len(urls)) + ' urls')
    return urls


def get_gaia_names():
    string = urlopen(GAIA_URL).read().decode('utf-8')
    pattern = re.compile('GaiaSource.*\.csv\.gz"')  # the pattern actually creates duplicates in the list
    filelist = pattern.findall(string)
    return list(map(lambda filename: filename[:-1], filelist))

def download(url):
    for i in range(0, 2):
        try:
            return pd.read_csv(url, compression='gzip', comment='#')
        except:
            logging.exception(f'Can\'t download ${url} (#{i})')
    return pd.read_csv(url, compression='gzip', comment='#')

def process(url, i):
    logging.info('Processing ' + str(i) + ': ' + url)
    df = download(url)
    df = df[df.radial_velocity.notnull()]
    df.to_csv('gdr3.csv', mode='a', header=i == 0)

def main():
    problem_urls_file = 'problem_urls.txt'

    urls = get_gaia_urls()
    for i in range(0, len(urls)):
        url = urls[i]
        try:
            process(url, i)
        except:
            logging.exception('Can\'t process dataset #' + str(i) + ': ' + url)
            with open(problem_urls_file, 'w') as f:
                f.write(str(i) + ' ' + url + '\n')
        print()



if __name__ == '__main__':
    main()
