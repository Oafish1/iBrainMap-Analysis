import os
import pickle

import pandas as pd


# Locations
DATA_FOLDER = '../../data/'
META = DATA_FOLDER + 'syn26527784_latest.csv'
CONTRAST = DATA_FOLDER + 'contrasts.csv'
DOSAGE = DATA_FOLDER + 'PsychAD_Dosage/genotype_varThresh0.5.csv'
COEX_FOLDER = DATA_FOLDER + 'freeze2/regulon_grn/'

# Freeze 2
# ATT = DATA_FOLDER + 'freeze2/attention/homo_5TF_1Tar_graph_with_edgeW_att.pkl'

# Freeze 2.5
ATT_FOLDER = DATA_FOLDER + 'freeze25/c01_5TF_10tar/'
ATT = ATT_FOLDER + 'output_embed_att/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_att.pkl'
GE = ATT_FOLDER + 'output_embed_att/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_graph_embedding.pkl'
SID = ATT_FOLDER + 'train_graph/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_sample_id.pkl'


### File Functions
def get_attention_columns():
    graphs_pkl = get_graphs_pkl()
    graph = graphs_pkl[list(graphs_pkl.keys())[0]]
    return [c for c in graph.columns if c not in ['from', 'to', 'from_gene', 'to_gene', 'att_mean', 'att_max']]


def load_graph_by_id(graph_id, source='attention', column=None, **kwargs):
    "Given a subject id `graph_id`, will return dataframe with column(s) `column` from `source`"
    # From individual graphs
    if source == 'coexpression':
        column = 'CoexWeight' if column is None else column
        # Get graph
        graph = pd.read_csv(f'{COEX_FOLDER}{graph_id}_regulon_list.csv')[['TF', 'gene', column, 'regulon']]
        graph = graph.rename(columns={'gene': 'TG', column: 'coef'})  # TF, TG, coef, regulon

    # From pkl
    elif source == 'attention':
        # columns
        # 'att_mean', 'att_max',
        # 'att_D_AD_0_1', 'att_D_AD_0_3', 'att_D_AD_0_5', 'att_D_AD_0_7',
        # 'att_D_no_prior_0', 'att_D_no_prior_1', 'att_D_no_prior_2', 'att_D_no_prior_3'
        column = 'att_max' if column is None else column  # Max for retention of head-specific prioritization

        # Load pkl
        graphs_pkl = get_graphs_pkl()

        # Get graph
        if type(column) == type([]):
            # If list, keep many columns and don't standardize name
            graph = graphs_pkl[graph_id][['from_gene', 'to_gene'] + column]
            graph = graph.rename(columns={'from_gene': 'TF', 'to_gene': 'TG'})
        else:
            graph = graphs_pkl[graph_id][['from_gene', 'to_gene', column]]
            graph = graph.rename(columns={'from_gene': 'TF', 'to_gene': 'TG', column: 'coef'})  # TF, TG, coef

    # Exception
    else:
        raise Exception(f'Source \'{source}\' not found.')

    return graph


def load_graph_embeddings():
    with open(SID, 'rb') as f:
        graph_sids = pickle.load(f)
    with open(GE, 'rb') as f:
        graph_embeddings = pickle.load(f)
    library = {sid: ge.detach().flatten().numpy() for sid, ge in zip(graph_sids, graph_embeddings)}
    return library


graphs_pkl = None
def get_graphs_pkl():
    # Load pkl if not already loaded
    global graphs_pkl
    if not graphs_pkl:
        with open(ATT, 'rb') as f:
            graphs_pkl = pickle.load(f)
    return graphs_pkl


def get_meta():
    return pd.read_csv(META)


def get_dosage():
    dosage = pd.read_csv(DOSAGE)
    dosage = dosage.set_index('snp_id')
    dosage = dosage.loc[dosage.index != 'snp_id']
    return dosage


contrast_table = None
def get_contrast(contrast):
    # Load if not already loaded
    global contrast_table
    if contrast_table is None:
        contrast_table = pd.read_csv(CONTRAST, dtype=str)

    # Construct dictionary
    library = {}
    for sub in pd.unique(contrast_table[contrast]):
        if pd.isna(sub): continue
        library[sub] = list(contrast_table['SubID'].loc[contrast_table[contrast]==sub])

    return library


def load_many_graphs(subject_ids, **kwargs):
    "Load as many graphs from `subject_ids` as available"
    graphs = []
    sids = []
    for sid in subject_ids:
        try:
            graphs.append(load_graph_by_id(sid, **kwargs))
            sids.append(sid)
        except: pass

    return graphs, sids


def get_enrichment(fname, num_diseases=10):
    # Skip if doesn't exist
    if not os.path.isfile(fname + '.csv'):
        return

    # Load file
    enrichment = pd.read_csv(fname + '.csv')
    enrichment = (
        enrichment[['Description', '_LogP_MyList']]
        .rename(columns={'Description': 'Disease', '_LogP_MyList': '-log(p)'})
    )
    enrichment['-log(p)'] *= -1
    enrichment['Cell Type'] = 'All'

    # Sort and filter to top 10
    # TODO: Filter to common diseases across cell types
    enrichment = enrichment.sort_values('-log(p)', ascending=False).iloc[:num_diseases]

    # Return
    return enrichment
