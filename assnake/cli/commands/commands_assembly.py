#
#
#

import click, os
@click.command('prepare_set_for_assembly')
@click.option('--df', '-d')
@click.option('--preproc','-p', help='Preprocessing to use', required = False)
@click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )
@click.option('--set_name','-n', help='Name of your sample set' )


def prepare_set_for_assembly(df, preproc, samples_to_add,set_name):
    print(df)
    print(preproc)
    print(samples_to_add)
    if preproc is None:
        preproc = 'raw'
    samples_to_add = [] if samples_to_add == '' else [c.strip() for c in samples_to_add.split(',')]
    df = assnake.api.loaders.load_df_from_db(df)
    ss = assnake.SampleSet.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)

    click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc', 'df']].sort_values('reads'), 
        headers='keys', tablefmt='fancy_grid'))

    set_dir = os.path.join(df['fs_prefix'], df['df'], 'assembly', set_name)
    os.makedirs(set_dir, exist_ok=True)

    set_loc = os.path.join(set_dir, 'sample_set.tsv')
    ss.samples_pd[['df', 'preproc', 'fs_name']].to_csv(set_loc, sep='\t', index=False)

