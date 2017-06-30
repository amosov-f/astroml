import errno
import os
import re
from urllib2 import urlopen

GAIA_URL = 'http://cdn.gea.esac.esa.int/Gaia/gaia_source/csv/'


def get_gaia_urls():
    return map(lambda filename: GAIA_URL + filename, get_gaia_names())


def get_gaia_names():
    string = urlopen(GAIA_URL).read()

    pattern = re.compile('GaiaSource.*\.csv\.gz"')  # the pattern actually creates duplicates in the list

    filelist = pattern.findall(string)

    return map(lambda filename: filename[:-1], filelist)


def get_gaia_filenames():
    return map(lambda url: get_gaia_filename(url), get_gaia_urls())


def get_gaia_filename(url):
    return os.sep.join(url.split(os.sep)[-3:])


def download_file(url):
    filename = get_gaia_filename(url)
    if os.path.exists(filename):
        print url + ' already exists in ' + filename
        return filename

    print 'Downloading ' + url
    try:
        make_dirs(filename)
        localfile = open(filename, 'wb')
        remotefile = urlopen(url)
        localfile.write(remotefile.read())
        localfile.close()
        remotefile.close()
        return filename
    except Exception, err:
        print "Can't download and save url " + url
        print err
        return None


def make_dirs(filename):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

                # response = urlopen(url)
                # if response.info().get('Content-Type') == 'application/x-gzip':
                #     buf = StringIO(response.read())
                #     f = gzip.GzipFile(fileobj=buf)
                #     data = f.read()
                # else:
                #     data = response.read()
                # response.close()
                # return data


                # urls = get_gaia_urls()
                # for url in urls:
                #     download_chunk(url)


                # for filename in filelist:
                #     filename = filename[:-1]
                #     file_url = url + filename
                #
                #     remotefile = urlopen(file_url)
                #     localfile = open('gaia_source/' + filename, 'wb')
                #     localfile.write(remotefile.read())
                #     localfile.close()
                #     remotefile.close()
