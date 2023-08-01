import colorsys
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
        alpha = get_alpha(g.ep.coef[e])
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
    """
    Subsets a graph to nodes in `g` connected to those in `vertices`.
    """
    # TODO: Add degrees of difference, i.e. 2nd gen connections and below
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


def concatenate_graphs(*graphs, color=True, exclude_zeroes_from_mean=True):
    """
    Concatenate all graphs provided, assesses duplicates by IDs
    """
    g = None
    g_coefs = {}
    for i, gc in enumerate(graphs):
        # Copy if first
        if not g:
            _add_attribute_to_dict(
                g_coefs,
                gc.edges(),
                indexer=lambda e: _get_edge_string(gc, e),
                attribute=gc.ep.coef,
                default=lambda: [0 for _ in range(len(graphs))],
                index=i,
            )
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


        # Track coefs
        _add_attribute_to_dict(
            g_coefs,
            gc.edges(),
            indexer=lambda e: _get_edge_string(gc, e),
            attribute=gc.ep.coef,
            default=lambda: [0 for _ in range(len(graphs))],
            index=i,
        )

        # Concatenate (assumes all vertex and edge properties are the same)
        # TODO: Other props, add coef averaging
        g, props = gt.graph_union(
            g,
            gc,
            intersection=vertex_map,
            props=[
                (g.vp.color, gc.vp.color),
                (g.vp.ids, gc.vp.ids),
                (g.vp.text, gc.vp.text),
                (g.vp.text_synthetic, gc.vp.text_synthetic),
                (g.ep.coef, gc.ep.coef),
            ],
            include=True)
        g.vp.color, g.vp.ids, g.vp.text, g.vp.text_synthetic, g.ep.coef = props

    # Label self loops
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop)

    # Add processed attributes
    g.vp.self_loop_value = g.new_vertex_property('double')
    g.ep.color = g.new_edge_property('vector<double>')
    for e in g.edges():
        # Get coefs
        coefs = g_coefs[_get_edge_string(g, e)]

        # Get processed attributes
        in_graph = [c != 0 for c in coefs]
        present_coef = np.mean([c for c in coefs if c != 0]) if exclude_zeroes_from_mean else np.mean(coefs)
        color = [0, 0, 0, 0]

        # Set color
        if sum(in_graph) == len(graphs) or not color:
            color[:3] = [0, 0, 0]
        else:
            # Get color index
            cindex = ''.join([str(int(b)) for b in in_graph])[::-1]  # Convert to binary
            cindex = int(cindex, 2) # Binary to int
            color[:3] = colorsys.hsv_to_rgb(cindex*(1./len(graphs)), 1, 1)

        # Set alpha
        color[3] = get_alpha(present_coef)

        # Write
        g.ep.color[e] = color
        if g.ep.self_loop:
            g.vp.self_loop_value[e.source()] = present_coef

    return g


def _get_edge_string(g, e):
    return f'{g.vp.ids[e.source()]}-{g.vp.ids[e.target()]}'


def _add_attribute_to_dict(dict, iterator, *, indexer=lambda x: x, attribute, default=lambda: [], index=None):
    """
    Add result `dict[indexer[i]] = attribute[i]` for `i` in `iterator` in-place
    """
    for item in iterator:
        if not indexer(item) in dict:
            dict[indexer(item)] = default()
        if index is None:
            dict[indexer(item)].append(attribute[item])
        else:
            dict[indexer(item)][index] = attribute[item]


def _normalize_dict_item_length(dict, length, default=0):
    """
    Normalize dictionary array length in-place
    """
    for k in dict.keys():
        if len(dict[k]) > length:
            raise IndexError('Dictionary entry too long.')
        elif len(dict[k]) < length:
            dict[k] += [default for _ in range(length - len(dict[k]))]


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
