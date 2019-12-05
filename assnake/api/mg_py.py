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