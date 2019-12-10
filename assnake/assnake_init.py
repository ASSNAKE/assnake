import click, sys, os, glob, yaml

@click.group()
@click.version_option()
@click.pass_context
def cli(ctx):
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = os.path.join(dir_of_this_file, '../snakemake/config.yml')
    config = {}
    with open(config_loc, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    wc_config_loc = os.path.join(dir_of_this_file, 'wc_config.yaml')
    wc_config = {}
    with open(wc_config_loc, 'r') as stream:
        try:
            wc_config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    ctx.obj = {'config': config, 'wc_config': wc_config}
    pass #Entry Point

def if __name__ == "__main__":
    cli()