import click
# from assnake.assnake_cli import pass_environment


@click.command('status', short_help='Shows smth.')
@click.pass_obj

def cli(config):

    if config.get('requests', None) is None:
        config['requests'] = ['1']
    else:
        config['requests'] += ['2']
    click.echo('We are inside status command!!!')