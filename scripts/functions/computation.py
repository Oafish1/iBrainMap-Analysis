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
def compute_statistics(meta, x, x_sub, filter=.8, **kwargs):
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
    # Detect synthetic vertices
    if hub_present:
        synthetic_vertices = list(graph[graph['TG'] == 'hub']['TF'])
    else:
        synthetic_vertices = [
            s for s in np.unique(graph['TG'])
            if (
                s != s.upper()
                or s in ['EN', 'IN', 'OPC', 'PC', 'VLMC', 'PVM', 'SMC']
                or sum([s.startswith(t) for t in ['EN_', 'IN_', 'CD8_']])
            )]
        synthetic_vertices += ['EN', 'IN', 'OPC', 'PC']
    cell_vertices = [v for v in synthetic_vertices if v != 'hub']

    # Filter to high-valued edges
    if filter:
        graph = graph[graph['coef'] >= graph['coef'].quantile(filter)]

    # Convert graph to graph_tools
    list_of_tuples = list(graph.itertuples(index=False, name=None))
    g = gt.Graph(list_of_tuples, hashed=True, eprops=[('coef', 'double')])

    # Label self loops and add color
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop)
    g.vp.self_loop_value = g.new_vertex_property('double')
    g.ep.color = g.new_edge_property('vector<double>')
    for e in g.edges():
        # Label self-loops
        if g.ep.self_loop:
            g.vp.self_loop_value[e.source()] = g.ep.coef[e]

        # Add color to edges
        x = g.ep.coef[e]
        x = np.log10(1+x)  # Log scaling
        # x = x**(1/3)  # Power scaling
        alpha = x / (1+x)
        alpha = .05 + .95 * alpha  # Add floor
        g.ep.color[e] = [0, 0, 0, alpha]

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
            is_tg = g_nosynthetic.get_in_degrees([v])[0] > 0
            if is_tf and not is_tg:
                g.vp.color[v] = rgba_to_hex(palette[2])
            elif not is_tf and is_tg:
                g.vp.color[v] = rgba_to_hex(palette[3])
            elif is_tf and is_tg:
                g.vp.color[v] = rgba_to_hex(palette[4])
            else:
                # Only connections from synthetic node
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


def subset_by_hub(g, vertices, verticies_are_ids=True, include_synthetic=False):
    # Convert ids to vertices
    vertices_new = []
    if verticies_are_ids:
        for vid in vertices:
            vertices_new.append(gt.find_vertex(g, g.vp.ids, vid)[0])
    vertices = vertices_new

    # Select vertices
    vfilt = []
    for v in vertices:
        vfilt.append(v)
        for vn in v.all_neighbors():
            if include_synthetic or not g.vp.text_synthetic[vn]:
                vfilt.append(vn)

    # Convert to proper format
    vfilt = [v in vfilt for v in g.vertices()]

    return gt.GraphView(
        g,
        vfilt=vfilt,
    )


def concatenate_graphs(*graphs):
    "Concatenate all graphs provided, assesses duplicates by IDs"
    g = None
    for gc in graphs:
        # Copy if first
        if not g:
            g = gc.copy()
            continue

        # Get common vertices
        vertex_map = gc.new_vertex_property('int')
        for v in gc.vertices():
            # Find corresponding v in g
            vid = gc.vp.ids[v]
            idx = gt.find_vertex(g, g.vp.ids, vid)  # Not really idx, actually list of vertex refs

            # Cases for finding
            if not idx:
                vertex_map[v] = -1
            elif len(idx) == 1:
                vertex_map[v] = int(idx[0])
            else:
                raise LookupError(f'ID \'{vid}\' has duplicate entries in \'g\'.')


        # Concatenate (assumes all vertex and edge properties are the same)
        # TODO: Other props
        g, props = gt.graph_union(
            g,
            gc,
            intersection=vertex_map,
            props=[(g.vp.ids, gc.vp.ids), (g.ep.coef, gc.ep.coef)],
            include=True)
        g.vp.ids, g.ep.coef = props

    return g


def get_graph_pos(g, scale=None):
    # Compute scale
    if not scale:
        scale = 20 * np.sqrt(g.num_vertices())

    # Compute layout
    # sfdp_layout(g), gt.arf_layout(g, max_iter=1000), radial_tree_layout(g, root), random_layout(g)
    # return gt.sfdp_layout(g, eweight=g.ep.coef)
    # return gt.arf_layout(g, weight=g.ep.coef)
    return gt.fruchterman_reingold_layout(g, weight=g.ep.coef, scale=scale)


def convert_vertex_map(source_graph, target_graph, vertex_map, debug=False):
    # NOTE: Probably a way to do this without `source_graph`
    converted_map = target_graph.new_vertex_property('vector<double>')
    for v in source_graph.vertices():
        # Find corresponding v in target_graph
        vid = source_graph.vp.ids[v]
        idx = gt.find_vertex(target_graph, target_graph.vp.ids, vid)

        # Cases for finding
        if idx:
            if len(idx) > 1: raise LookupError(f'ID \'{vid}\' has duplicate entries in \'g\'.')
            converted_map[idx[0]] = vertex_map[v]

    if debug:
        for v in target_graph.vertices():
            if not converted_map[v]:
                print(f'\'{target_graph.vp.ids[v]}\' not defined in `vertex_map`.')

    return converted_map
