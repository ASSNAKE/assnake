import random
import string

def random_file_name():
    length = random.randint(2, 10)
    return ''.join([random.choice(string.ascii_lowercase) for _ in range(length)])


def random_path():
    nesting = random.randint(2, 10)
    return '/'.join([random_file_name() for _ in range(nesting)])