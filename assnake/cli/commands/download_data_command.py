import click
from assnake.utils import pathizer, download_from_url, get_url_remote_file_size, read_yaml
import os
import urllib

# ---------------------------------------------------------------------------------------
#                                   Download Data
# ---------------------------------------------------------------------------------------
@click.command(name='download-data')
@click.option('--output-dir', '-o', required = False, help="Output directory")
@click.option('--data','-d', prompt="Type in the data name to download",help="The data to download", type=click.Choice(['centrifuge', 'test', 'custom']))
@click.option('--force', '-f', is_flag=True, help='Do it without questions or confirmations')
@click.pass_obj
def data_download(context, output_dir, data,  force):
    """
    Download data from widely-used predefined databases or from your custom url
    """
    dir_of_file = os.path.dirname(os.path.abspath(__file__))
    config = read_yaml(os.path.join(dir_of_file, 'config_download_data_command.yaml'))
    path = output_dir
    if path is None and not force:
        if click.confirm('Current directory will be used. Continue?', abort=True):
            path = ''
    path = pathizer(path)
    try:
        os.makedirs(path)
    except FileExistsError as e:
        if not force:
            click.echo("Existing directory will be used")
    else:
        if not force:
            click.secho("New directory "+click.style(path, fg='blue')+ " was created")

    if data == 'custom':
        url = click.prompt('Type in yout url')
    else:

        url = config[data]
    try:
        size = get_url_remote_file_size(url)
    except urllib.error.HTTPError as e:
        click.secho("HTTP error"+click.style(e.code, fg='red'))
        exit(2)
    if force or click.confirm("The {} bytes will be doawnloaded, continue?".format(size), abort=True):
        try:
            size = download_from_url(url, os.path.join(path, config["{}_name".format(data)]))
        except urllib.error.HTTPError as e:
            click.secho("HTTP error" + click.style(e.code, fg='red'))
            exit(2)
    click.secho("Successfully downloaded", bg='green', fg='bright_white')
    click.echo("Directory {}, size {}".format(path, size))