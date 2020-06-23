import yaml, configparser, os, click
import os, sys
import requests, urllib
from tqdm import tqdm

from pathlib import Path

def read_yaml(file_location):
    yaml_file = {}
    with open(file_location, 'r') as stream:
        try:
            yaml_file = yaml.load(stream, Loader=yaml.FullLoader)
            return yaml_file
        except yaml.YAMLError as exc:
            print(exc)

def load_wc_config():
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    return (read_yaml(os.path.join(dir_of_this_file, './snake/wc_config.yaml')))


def read_internal_config():
    '''
    Reads the config at ~/.config/assnake/internal_config.yaml and returns as dict
    '''
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    
    if os.path.isfile(internal_config_loc):
        internal_config = read_yaml(internal_config_loc)
        return internal_config
    else:
        internal_config_dir = os.path.join(str(Path.home()), '.config/assnake/')

        os.makedirs(internal_config_dir, exist_ok=True)

        if not os.path.isfile(internal_config_loc):
            with open(internal_config_loc, 'w+') as file:
                _ = yaml.dump({'assnake_db' : '', 'instance_config_loc': ''}, file, sort_keys=False)
        return {'assnake_db' : '', 'instance_config_loc': ''}

def read_assnake_instance_config():
    '''
    Reads particular assnake instance config. It is stored inside assnake database as config.yaml (Name subject to change). 
    :return: Returns dict if instance config exists, None otherwise.
    '''
    internal_config = read_internal_config()
    instance_config_loc = os.path.join(internal_config['assnake_db'], 'config.yaml')

    if os.path.isfile(instance_config_loc):
        return read_yaml(instance_config_loc)
    else:
        return None

def update_internal_config(update_dict):
    '''
    Updates internal config with provided dictionary

    :return: updated internal config if everything is ok, None otherwise
    '''
    internal_config = read_internal_config()
    internal_config.update(update_dict)
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    with open(internal_config_loc, 'w+') as file:
        _ = yaml.dump(internal_config, file, sort_keys=False)

    return internal_config

def check_if_assnake_is_initialized():
    '''
    Just tries to read assnake instance config, and prints an error if there is no instance config.
    '''
    instance_config = read_assnake_instance_config()

    if instance_config is None:
        click.secho("You need to init your installation!", fg='red', bold=True)
        click.echo("Don't worry, it won't take long.")
        click.echo('Just run ' + click.style('assnake init start', bg='blue', fg='bright_white'))
        exit()




## ====RENAME THIS === update_instance_config
def update_config(update_dict):
    '''
    Updates instance config with provided dictionary

    :return: updated instance config if everything is ok, None otherwise
    '''

    instance_config = read_assnake_instance_config()
    instance_config.update(update_dict)

    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')

    with open(internal_config_loc, 'w+') as file:
        _ = yaml.dump(instance_config, file, sort_keys=False)

    return instance_config


## ============================================






def pathizer(path):
    """
    Make from path absolute path
    :param path: absolute or relative path
    :return: absolute path
    """
    if path is None or path==' ' or path == '':
        return os.popen('pwd').read().replace('\n', '')
    final = '{prefix}/{rel_path}'.format(prefix=os.popen('pwd').read().replace('\n', ''), rel_path=path) if (path[0] != '/') else path
    return final.rstrip('\/')


def dict_norm_print(d, indent=1):
    if isinstance(d, str):
        click.echo('\t' * indent + str(d), nl=False)
    else:
        for key, value in d.items():
            click.echo('\t' * indent + str(key), nl=False)
            if isinstance(value, dict):
                click.echo('--↴	')
                dict_norm_print(value, indent + 1)
            else:
                click.echo('\t' * (indent) + str(value))


## {{{ http://code.activestate.com/recipes/578019/ (r15)
# see: http://goo.gl/kTQMs
SYMBOLS = {
    'customary': ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
    'customary_ext': ('byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa',
                      'zetta', 'iotta'),
    'iec': ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
    'iec_ext': ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                'zebi', 'yobi'),
}


def bytes2human(n, format='%(value).1f %(symbol)s', symbols='customary'):
    """
    Convert n bytes into a human readable string based on format.
    symbols can be either "customary", "customary_ext", "iec" or "iec_ext",
    see: http://goo.gl/kTQMs

      >>> bytes2human(0)
      '0.0 B'
      >>> bytes2human(0.9)
      '0.0 B'
      >>> bytes2human(1)
      '1.0 B'
      >>> bytes2human(1.9)
      '1.0 B'
      >>> bytes2human(1024)
      '1.0 K'
      >>> bytes2human(1048576)
      '1.0 M'
      >>> bytes2human(1099511627776127398123789121)
      '909.5 Y'

      >>> bytes2human(9856, symbols="customary")
      '9.6 K'
      >>> bytes2human(9856, symbols="customary_ext")
      '9.6 kilo'
      >>> bytes2human(9856, symbols="iec")
      '9.6 Ki'
      >>> bytes2human(9856, symbols="iec_ext")
      '9.6 kibi'

      >>> bytes2human(10000, "%(value).1f %(symbol)s/sec")
      '9.8 K/sec'

      >>> # precision can be adjusted by playing with %f operator
      >>> bytes2human(10000, format="%(value).5f %(symbol)s")
      '9.76562 K'
    """
    n = int(n)
    if n < 0:
        raise ValueError("n < 0")
    symbols = SYMBOLS[symbols]
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i + 1) * 10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    return format % dict(symbol=symbols[0], value=n)


def human2bytes(s):
    """
    Attempts to guess the string format based on default symbols
    set and return the corresponding bytes as an integer.
    When unable to recognize the format ValueError is raised.

      >>> human2bytes('0 B')
      0
      >>> human2bytes('1 K')
      1024
      >>> human2bytes('1 M')
      1048576
      >>> human2bytes('1 Gi')
      1073741824
      >>> human2bytes('1 tera')
      1099511627776

      >>> human2bytes('0.5kilo')
      512
      >>> human2bytes('0.1  byte')
      0
      >>> human2bytes('1 k')  # k is an alias for K
      1024
      >>> human2bytes('12 foo')
      Traceback (most recent call last):
          ...
      ValueError: can't interpret '12 foo'
    """
    init = s
    num = ""
    while s and s[0:1].isdigit() or s[0:1] == '.':
        num += s[0]
        s = s[1:]
    num = float(num)
    letter = s.strip()
    for name, sset in SYMBOLS.items():
        if letter in sset:
            break
    else:
        if letter == 'k':
            # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
            sset = SYMBOLS['customary']
            letter = letter.upper()
        else:
            raise ValueError("can't interpret %r" % init)
    prefix = {sset[0]: 1}
    for i, s in enumerate(sset[1:]):
        prefix[s] = 1 << (i + 1) * 10
    return int(num * prefix[letter])

## end of http://code.activestate.com/recipes/578019/ }}}


## https://gist.github.com/wy193777/0e2a4932e81afc6aa4c8f7a2984f34e2
def download_from_url(url, dst):
    """
    @param: url to download file
    @param: dst place to put the file
    """

    req = urllib.request.Request(url,  method='HEAD')
    f = urllib.request.urlopen(req)

    file_size = int(f.headers['Content-Length'])
    if os.path.exists(dst):
        first_byte = os.path.getsize(dst)
    else:
        first_byte = 0
    if first_byte >= file_size:
        return file_size
    header = {"Range": "bytes=%s-%s" % (first_byte, file_size)}
    pbar = tqdm(
        total=file_size, initial=first_byte,
        unit='B', unit_scale=True, desc=url.split('/')[-1])
    req = requests.get(url, headers=header, stream=True)
    with(open(dst, 'ab')) as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                pbar.update(1024)
    pbar.close()
    return file_size