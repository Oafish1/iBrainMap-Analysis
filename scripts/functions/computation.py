import itertools as it

import graph_tool.all as gt
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from .file import *
from .utility import *


### Computational Functions
def compute_statistics(meta, x, x_sub, filter=.8, **kwargs):
    # Initialize graph summary list
    graph_summary_list = []

    # Get unique values
    unique_x = meta[x].unique()
    unique_x_sub = meta[x_sub].unique()

    # Calculate per graph
    print('Calculating statistics...')
    for val, val_sub in tqdm(it.product(unique_x, unique_x_sub), total=len(unique_x)*len(unique_x_sub)):
        graph_ids = list(meta[(meta[x]==val)*(meta[x_sub]==val_sub)]['SubID'])

        num_samples = len(graph_ids)  # 10
        for graph_id in np.random.choice(graph_ids, min(len(graph_ids), num_samples), replace=False):
            # Load graph
            try:
                graph = load_graph_by_id(graph_id, **kwargs)
            except:
                continue

            # Filter synthetic cells
            if 'source' in kwargs and kwargs['source'] == 'attention':
                synthetic_vertices = list(graph[graph['TG'] == 'hub']['TF'])
            else:
                synthetic_vertices = []

            # TODO: Filter by coef, cell-type regulons

            # Results
            results_list = [graph_id, val, val_sub]

            # Filter to high coef
            if filter:
                graph = graph[graph['coef'] >= graph['coef'].quantile(filter)]

            # Calculate TF outgoing
            tf_outgoing = np.mean([(graph['TF']==gene).sum() for gene in graph['TF'].unique() if gene not in synthetic_vertices])
            results_list.append(tf_outgoing)

            # Calculate TG outgoing (Long Runtime)
            tg_outgoing = 0  # np.mean([(graph['TG']==gene).sum() for gene in graph['TG'].unique() if gene not in synthetic_vertices])
            results_list.append(tg_outgoing)

            # Create graph
            graph_nx = nx.from_pandas_edgelist(graph, 'TF', 'TG', 'coef')

            # Calculate TF closeness
            tf_closeness = np.mean([nx.closeness_centrality(graph_nx, u=gene) for gene in graph['TF'].unique() if gene not in synthetic_vertices])
            results_list.append(tf_closeness)

            # Get cliques
            cliques = sum([1 for _ in nx.find_cliques(graph_nx)])
            results_list.append(cliques)

            # Record results
            graph_summary_list.append(results_list)

    # Create summary df
    graph_summary = pd.DataFrame(graph_summary_list, columns=[
        'Graph ID',
        x,
        x_sub,
        'TF Outgoing',
        'TG Outgoing',
        'TF Closeness',
        'Cliques',
    ])

    # Return
    return graph_summary


def compute_graph(graph, filter=0, hub_present=False):  # TODO: Find something more interesting than cutting off `filter`
    # Filter to high-valued edges
    if filter:
        graph = graph[graph['coef'] >= graph['coef'].quantile(filter)]

    # Convert graph to graph_tools
    list_of_tuples = list(graph.itertuples(index=False, name=None))
    g = gt.Graph(list_of_tuples, hashed=True, eprops=[('coef', 'double')])

    # Label self loops and add color
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop, mark_only=True)
    g.vp.self_loop_value = g.new_vertex_property('double')
    g.ep.color = g.new_edge_property('vector<double>')
    for e in g.edges():
        # Label self-loops
        if g.ep.self_loop[e]:
            g.vp.self_loop_value[e.source()] = g.ep.coef[e]

        # Add color to edges
        alpha = get_alpha(g.ep.coef[e])
        g.ep.color[e] = [0, 0, 0, alpha]

    # Determine color and flavor text
    g = assign_vertex_properties(g)

    # View without self-loops
    g_noself = gt.GraphView(g, efilt=lambda e: not g.ep.self_loop[e])

    return g_noself


def compute_edge_summary(graphs=None, concatenated_graph=None, *, subject_ids, min_common_edges=1):
    # Setup
    assert graphs is not None or concatenated_graph is not None
    if graphs is not None: concatenated_graph = concatenate_graphs(*graphs)

    # Format edge weights
    # df = pd.DataFrame(columns=['Edge']+subject_ids)
    df = {k: [] for k in ['Edge']+subject_ids}
    print('Collecting edges...')
    for e in tqdm(concatenated_graph.edges(), total=concatenated_graph.num_edges()):
        edge_name = get_edge_string(concatenated_graph, e)
        coefs = concatenated_graph.ep.coefs[e]
        # Take only edges which are common between two or more graphs
        if sum([c!=0 for c in coefs]) >= min_common_edges:
            row = [edge_name] + list(coefs)
            # df.loc[df.shape[0]] = row  # Slow
            for k, v in zip(df, row):
                df[k].append(v)
    df = pd.DataFrame(df)

    # Find variance and mean
    df['Variance'] = np.var(df.iloc[:, 1:1+len(coefs)], axis=1)
    df['Mean'] = np.mean(df.iloc[:, 1:1+len(coefs)], axis=1)
    # df = df.sort_values(['Variance', 'Mean'], ascending=[True, False])
    # df['index'] = list(range(len(df)))  # Store sort

    return df, concatenated_graph


def compute_aggregate_edge_summary(contrast_subject_ids, *, column, max_graphs=np.inf):
    # For each subgroup of the contrast
    contrast_concatenated_graphs = {}; contrast_concatenated_subject_ids = {}
    for key, subject_ids in contrast_subject_ids.items():  # , total=len(contrast_subject_ids)
        # Get concatenated graph
        graphs = []; sids = []
        for sid in np.random.choice(subject_ids, len(subject_ids), replace=False):  # , total=len(subject_ids)
            try: graphs.append(compute_graph(load_graph_by_id(sid, column=column)))
            except: continue
            sids.append(sid)
            if len(sids) >= max_graphs: break
        concatenated_graph = concatenate_graphs(*graphs)
        contrast_concatenated_graphs[key] = concatenated_graph
        contrast_concatenated_subject_ids[key] = sids

        # Cleanup
        del graphs

    return contrast_concatenated_graphs, contrast_concatenated_subject_ids
