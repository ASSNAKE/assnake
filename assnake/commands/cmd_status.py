import click
# from assnake.assnake_cli import pass_environment


@click.command('status', short_help='Shows smth.')
def cli():
    click.echo('Dynamically loaded command at runtime')