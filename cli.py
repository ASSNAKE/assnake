import click, sys
print(sys.path)

import api.loaders

@click.group()
def datasets():
    pass

@click.command()
def df_list():
    print(api.loaders.load_dfs_from_db(''))

datasets.add_command(df_list)

if __name__ == '__main__':
    datasets()