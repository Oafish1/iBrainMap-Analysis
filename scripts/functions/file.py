import pickle

import pandas as pd


DATA_FOLDER = '../../data/'
META = DATA_FOLDER + 'syn26527784_latest.csv'
CONTRAST = DATA_FOLDER + 'contrasts.csv'
COEX_FOLDER = DATA_FOLDER + 'freeze2/regulon_grn/'
FREEZE2_ATT = DATA_FOLDER + 'freeze2/attention/homo_5TF_1Tar_graph_with_edgeW_att.pkl'
FREEZE25_ATT_FOLDER = DATA_FOLDER + 'freeze25/c01_5TF_10tar/'
FREEZE25_ATT = FREEZE25_ATT_FOLDER + 'output_embed_att/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_att.pkl'
FREEZE25_GE = FREEZE25_ATT_FOLDER + 'output_embed_att/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_graph_embedding.pkl'
FREEZE25_SID = FREEZE25_ATT_FOLDER + 'train_graph/homo_c01_5TF_10Tar_p2_5_graph_with_edgeW_sample_id.pkl'


### File Functions
graphs_pkl = None
def load_graph_by_id(graph_id, source='attention', column=None, **kwargs):
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
        column = 'att_mean' if column is None else column

        # Load pkl if not already loaded
        global graphs_pkl
        if not graphs_pkl:
            with open(FREEZE25_ATT, 'rb') as f:
                graphs_pkl = pickle.load(f)
        # Get graph
        graph = graphs_pkl[graph_id][['from_gene', 'to_gene', column]]
        graph = graph.rename(columns={'from_gene': 'TF', 'to_gene': 'TG', column: 'coef'})  # TF, TG, coef

    # Exception
    else:
        raise Exception(f'Source \'{source}\' not found.')

    return graph


def load_graph_embeddings():
    with open(FREEZE25_SID, 'rb') as f:
        graph_sids = pickle.load(f)
    with open(FREEZE25_GE, 'rb') as f:
        graph_embeddings = pickle.load(f)
    library = {sid: ge.detach().flatten().numpy() for sid, ge in zip(graph_sids, graph_embeddings)}
    return library


def get_meta():
    return pd.read_csv(META)


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
