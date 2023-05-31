import errno
import gzip
import os
import re
from io import StringIO
from typing import Optional
from urllib.request import urlopen
import pandas as pd

GAIA_URL = 'http://cdn.gea.esac.esa.int/Gaia/gdr2/vari_cepheid/csv/'


def get_gaia_urls():
    return map(lambda filename: GAIA_URL + filename, get_gaia_names())


def get_gaia_names():
    string = urlopen(GAIA_URL).read().decode('utf-8')

    pattern = re.compile('VariCepheid.*\.csv\.gz"')  # the pattern actually creates duplicates in the list

    filelist = pattern.findall(string)

    return map(lambda filename: filename[:-1], filelist)


def get_gaia_filenames():
    return map(lambda url: get_gaia_filename(url), get_gaia_urls())


def get_gaia_filename(url):
    return os.sep.join(url.split(os.sep)[-3:])


def download_file(url) -> Optional[str]:
    filename = get_gaia_filename(url)
    if os.path.exists(filename):
        print(url + ' already exists in ' + filename)
        return filename

    print('Downloading ' + url)
    try:
        make_dirs(filename)
        localfile = open(filename, 'wb')
        remotefile = urlopen(url)
        localfile.write(remotefile.read())
        localfile.close()
        remotefile.close()
        return filename
    except:
        print("Can't download and save url " + url)
        return None


def make_dirs(filename):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def main():
    result = None
    for url in get_gaia_urls():
        filename = download_file(url)
        with gzip.open(filename, 'rb') as f:
            file_content = f.read().decode('utf-8')
            df = pd.read_csv(StringIO(file_content))
            if result is None:
                result = df
            else:
                result = pd.concat((result, df))
    print(result)
    result.to_csv('VariCepheid.csv')



if __name__ == '__main__':
    main()