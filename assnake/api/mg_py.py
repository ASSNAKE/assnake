import pandas as pd

def create_taxa_count_from_dada2(tax_table, otu_table, rank_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']):
    tax_table = tax_table.copy().sort_index(inplace = False, axis=0)
    otu_table = otu_table.copy().sort_index(inplace = False, axis=0)

    cnt = otu_table.sort_index(inplace = False).copy()
    # sort columns
    cnt = cnt.reindex(sorted(cnt.columns), axis=1)

    # for multi-index
    count_nd_taxa = pd.concat([tax_table, cnt], axis=1)
    # filling NAs with underscore, for better names construction
    # count_nd_taxa = count_nd_taxa.fillna('_')
    count_nd_taxa[rank_names] = count_nd_taxa[rank_names].fillna('_')
    # available rank_names
    rank_names = list(tax_table.columns)
    # Actual dataframe
    taxa_counts = count_nd_taxa.set_index(rank_names, inplace=False)

    return taxa_counts

def tax_glom(taxa_counts, rank='Phylum', include_na = True, index = 'simple_long', 
    rank_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
    ranks_short = ['k', 'p', 'c', 'o', 'f', 'g', 's']):
    '''
    Agglomerates data at desired taxonomic rank
    
    Args:
        tax_table (:obj:`pandas.DataFrame`): DataFrame with counts and multiindex with taxonomic information. Samples are columns. 
        rank (str): Rank level at wich we want to agglomerate data
        include_na (bool): Whether to include counts without classification at desired level
    
    Returns:
        agg_counts (:obj:`pandas.DataFrame`):
    '''
    
    # rank_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    # ranks_short = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    
    ranks = rank_names[:rank_names.index(rank)+1]
    print(ranks)
    agg_counts = taxa_counts.copy()
    agg_counts = taxa_counts.groupby(by = ranks).sum()
    agg_counts = agg_counts.reset_index() 
    
    if index == 'simple_long':
        # ranks_short = ['k', 'p', 'c', 'o', 'f', 'g', 's']
        for i, t in enumerate(ranks):
            if i == 0:
                agg_counts['OTU'] = ranks_short[i]+'__' + agg_counts[t]
            else:
                agg_counts['OTU'] = agg_counts['OTU'] + ';'+ranks_short[i]+'__'+agg_counts[t]
                
        agg_counts.index = agg_counts['OTU']
        agg_counts = agg_counts.drop(ranks + ['OTU'], axis=1)

    elif index == 'simple_short':
        print('simple_short')
    elif index == 'multiindex':
        agg_counts = agg_counts.set_index(ranks, inplace=False)
    
    return agg_counts

def single_level_from_mp2_table(tax_table, rank = 'g__', ranks = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 't__']):
    """
    Takes metaphlan2 style table with all taxonomic levels included, like
    k__Archaea; k__Archaea|p__Euryarchaeota; etc., and desired taxonomic level as input, returns feature table 
    with only selected taxonomic level features.

    Args:
        tax_table (:obj:`pandas.DataFrame`): DataFrame with counts and full taxonomic information. Samples are columns. 
        rank (str): Rank level at wich we want to agglomerate data
        levels (list(str)): List with all ranks present in table
    Returns:
        tax_table_pruned (:obj:`pandas.DataFrame`): 
    """

    ind = ranks.index(rank)
    columns = tax_table.columns
    cols = []
    for col in columns:
        if not (ranks[ind] in col and ranks[ind+1] not in col):
            cols.append(col)
    if rank == 'k__':
        cols.remove('UNKNOWN')
    mp2_all_order = tax_table.drop(cols, axis=1)
    mp2_all_order = mp2_all_order.fillna(0)

    mm = mp2_all_order.div(mp2_all_order.sum(axis=1), axis=0)
    mm = mm.fillna(0).T

    return mm