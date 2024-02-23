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

def print_python_help(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return

    command_name = ctx.command.name
    excluded_params = ['custom_help', 'show_python', 'samples_to_add', 'exclude_samples', 'preprocessing']

    python_snippet = f"# Python configuration for {command_name}\n"
    python_snippet += f"{command_name}_params = {{\n"

    for option in ctx.command.params:
        if isinstance(option, click.Option) and option.name not in excluded_params:
            param_value = ctx.params.get(option.name, 'your_value_here')
            # Add a comment only if the default value is used
            comment = '' if option.name in ctx.params else ''
            python_snippet += f"    '{option.name}': '{param_value}',{comment}\n"

    python_snippet += "}\n"
    python_snippet += f"pipeline.add_result('{command_name}', **{command_name}_params)\n"

    click.echo(python_snippet)
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
              help='Show this message and exit.'),
    click.option('--show-python', 
                 is_flag=True, callback=print_python_help, 
                 help=("Generate a Python snippet for configuring this command in a Python pipeline. "))

]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options

