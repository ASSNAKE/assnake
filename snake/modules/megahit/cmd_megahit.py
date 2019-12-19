import click, glob, os
import assnake.api.loaders
import assnake.api.sample_set


@click.group('megahit', short_help='Ultra-fast and memory-efficient NGS assembler')
def cli():
    pass

@click.command('mh', short_help='Start the assembler')
@click.option('--df','-d', help='Name of the dataset', required=True )
@click.pass_obj

def mh(config, df):
    res_list = []
    df = assnake.api.loaders.load_df_from_db(df)
    prepared_sets = glob.glob(os.path.join(df['fs_prefix'],df['df'],'assembly/*/*/sample_set.tsv'))
    click.echo('We found ' + str(len(prepared_sets)) + ' prepared sets:')
    for i, pset in enumerate(prepared_sets):
        name = pset.split('/')[-2]
        run_info = pset.split('/')[-3]
        click.echo(click.style(str(i), bold = True) + ' ' + run_info + ' ' + name)

    selected_sets = []
    sel_sets_str = click.prompt('Please, enter sets you want to assemble, comma separated. You can also type all, to select all sets', type=str)
    sel_sets_str = sel_sets_str.replace(' ', '')
    if sel_sets_str != 'all':
        selected_sets = [int(s) for s in sel_sets_str.split(',')] 
    click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
    min_len = click.prompt('Please, enter minimum contig lenth', type=int)

    for ss in selected_sets:
        res_list += [prepared_sets[ss].replace('sample_set.tsv', 'final_contigs__{min_len}.fa'.format(min_len=min_len))]

cli.add_command(mh)