import click
# from assnake.assnake_cli import pass_environment


@click.command('sample-set-for-assembly', short_help='Select what samples you want to use for (co)assembly')
@click.option('--df','-d', help='Main dataset' )
@click.pass_obj
def cli(config, df):
    print(config)
    click.echo(df)