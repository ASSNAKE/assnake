import click

# https://stackoverflow.com/a/40195800
sample_set_construction_options = []

sample_set_construction_options = [
    click.option('--df','-d', help='Name of the dataset'),
    click.option('--preproc','-p', help='Preprocessing to use' ),

    click.option('--meta-column', '-c', help='Select samples based on metadata column' ),
    click.option('--column-value','-v', help='Value of metadata column by which select samples' ),

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
]

options_w_params = [
    click.option('--df','-d', help='Name of the dataset', required=True ),
    click.option('--preproc','-p', help='Preprocessing to use' ),

    click.option('--meta-column', '-c', help='Select samples based on metadata column' ),
    click.option('--column-value','-v', help='Value of metadata column by which select samples' ),

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
    click.option('--params', 
                help='Parameters to use', 
                default='def',
                type=click.STRING )
]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options

