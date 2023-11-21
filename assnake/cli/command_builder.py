import click


def custom_help(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    dataset = ctx.params.get('df')
    print(ctx.params)
    if dataset:
        print('dataset provided')
        preset_option = next((option for option in ctx.command.params if option.name == 'preset'), None)
        if preset_option:
            print('We have preset option')

            preset_option.help = f"Choose a preset for the command. Available presets for dataset {dataset}: [preset1, preset2]"

    # Generate the help message based on the dataset
    help_message = f"Help message for dataset: {dataset}"
    click.echo(help_message)

    click.echo(ctx.get_help())

    ctx.exit()

# https://stackoverflow.com/a/40195800
sample_set_construction_options = []

sample_set_construction_options = [
    click.option('--dataset','-d', help='Name of the dataset'),
    click.option('--preprocessing','-p', help='Preprocessing to use' ),

    # click.option('--meta-column', '-c', help='Select samples based on metadata column' ),
    # click.option('--column-value','-v', help='Value of metadata column by which select samples' ),

    click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
    click.option('--exclude-samples','-x', 
                help='Exclude this samples from run', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
    click.option('--custom-help', is_flag=True, callback=custom_help, 
              expose_value=False, 
              help='Show this message and exit.')
]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options

