import pickle

import pandas as pd


### File Functions
graphs_pkl = None
def load_graph_by_id(graph_id, source='attention', column=None, **kwargs):
    # From individual graphs
    if source == 'coexpression':
        column = 'CoexWeight' if column is None else column
        # Get graph
        graph = pd.read_csv(f'../../data/PsychAD_freeze2_personalized_grpahs/regulon_grn/{graph_id}_regulon_list.csv')[['TF', 'gene', column, 'regulon']]
        graph = graph.rename(columns={'gene': 'TG', column: 'coef'})  # TF, TG, coef, regulon

    # From pkl
    elif source == 'attention':
        column = 'att_mean' if column is None else column
        # columns
        # 'att_mean', 'att_max',
        # 'att_D_AD_0_1', 'att_D_AD_0_3', 'att_D_AD_0_5', 'att_D_AD_0_7',
        # 'att_D_no_prior_0', 'att_D_no_prior_1', 'att_D_no_prior_2', 'att_D_no_prior_3'
        # Load pkl if not already loaded
        global graphs_pkl
        if not graphs_pkl:
            with open(f'../../data/ting/2023-06-26/homo_5TF_1Tar_graph_with_edgeW_att.pkl', 'rb') as f:
                graphs_pkl = pickle.load(f)
        # Get graph
        graph = graphs_pkl[graph_id][['from_gene', 'to_gene', column]]
        graph = graph.rename(columns={'from_gene': 'TF', 'to_gene': 'TG', column: 'coef'})  # TF, TG, coef

    # Exception
    else:
        raise Exception(f'Source \'{source}\' not found.')

    return graph
