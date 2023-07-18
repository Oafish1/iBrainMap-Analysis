import itertools as it

import graph_tool.all as gt
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from .file import *
from .utility import *


### Computational Functions
def compute_statistics(meta, x, x_sub, filter=True, **kwargs):
    # Initialize graph summary list
    graph_summary_list = []

    # Get unique values
    unique_x = meta[x].unique()
    unique_x_sub = meta[x_sub].unique()

    # Calculate per graph
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
                graph = graph[graph['coef'] > graph['coef'].quantile(.8)]

            # Calculate TF outgoing
            tf_outgoing = np.mean([(graph['TF']==gene).sum() for gene in graph['TF'].unique() if gene not in synthetic_vertices])
            results_list.append(tf_outgoing)

            # Calculate TG outgoing (Long Runtime)
            tg_outgoing = 0  # np.mean([(graph['TG']==gene).sum() for gene in graph['TG'].unique() if gene not in synthetic_vertices])
            results_list.append(tg_outgoing)

            # Calculate TF closeness
            graph_nx = nx.from_pandas_edgelist(graph, 'TF', 'TG', 'coef')
            tf_closeness = np.mean([nx.closeness_centrality(graph_nx, u=gene) for gene in graph['TF'].unique() if gene not in synthetic_vertices])
            results_list.append(tf_closeness)

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
    ])

    # Return
    return graph_summary


def compute_graph(graph):
    # Detect synthetic vertices
    synthetic_vertices = list(graph[graph['TG'] == 'hub']['TF'])
    cell_vertices = [v for v in synthetic_vertices if v != 'hub']

    # Convert graph to graph_tools
    list_of_tuples = list(graph.itertuples(index=False, name=None))
    g = gt.Graph(list_of_tuples, hashed=True, eprops=[('coef', 'double')])

    # Label self loops and add vertex property
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop)
    g.vp.self_loop_value = g.new_vertex_property('double')
    for e in g.edges():
        if g.ep.self_loop:
            g.vp.self_loop_value[e.source()] = g.ep.coef[e]

    # View without synthetic nodes or self loops
    # Need to do `vfilt` slowly bc `g.vp.ids.fa`` doesn't work with string
    g_nosynthetic = gt.GraphView(
        g,
        vfilt=[g.vp.ids[v] not in synthetic_vertices for v in g.vertices()],
        efilt=1-g.ep.self_loop.fa,
    )

    # Determine color and flavor text
    # g.vp.color = g.new_vertex_property('vector<double>')
    g.vp.color = g.new_vertex_property('string')  # Can't show with text if set to `vector<double>``
    g.vp.text_synthetic = g.new_vertex_property('string')
    g.vp.text = g.new_vertex_property('string')
    for v in g.vertices():
        v_id = g.vp.ids[v]
        palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # Hub
        if v_id in ['hub']:
            g.vp.color[v] = rgba_to_hex(palette[0])
            g.vp.text_synthetic[v] = v_id
            g.vp.text[v] = v_id
            # root = v
        # Cell-type
        elif v_id in cell_vertices:
            g.vp.color[v] = rgba_to_hex(palette[1])
            g.vp.text_synthetic[v] = v_id
            g.vp.text[v] = v_id
        # Default
        else:
            is_tf = g_nosynthetic.get_out_degrees([v])[0] > 0
            is_tg = g_nosynthetic.get_in_degrees([v]) > 0
            if is_tf and not is_tg:
                g.vp.color[v] = rgba_to_hex(palette[2])
            elif not is_tf and is_tg:
                g.vp.color[v] = rgba_to_hex(palette[3])
            elif is_tf and is_tg:
                g.vp.color[v] = rgba_to_hex(palette[4])
            else:
                # raise Exception('How?')
                g.vp.color[v] = '#FFFFFF'
            # Add text to outer nodes (optional)
            g.vp.text[v] = v_id

    # View without self-loops
    g_noself = gt.GraphView(g, efilt=1-g.ep.self_loop.fa)

    return g_noself


def subset_graph(source, target):
    "Subset `source` graph nodes to those present in `target`"
    return gt.GraphView(
        source,
        vfilt=[source.vp.ids[v] in [target.vp.ids[v]
                                    for v in target.vertices()]
                for v in source.vertices()],
    )
