import colorsys

import graph_tool.all as gt
import numpy as np

from .file import *


### Utility functions
def rgba_to_hex(rgba):
    int_rgba = [int(255*i) for i in rgba]
    return '#{:02x}{:02x}{:02x}'.format(*int_rgba)


def rgba_array_to_rgba_string(array):
    int_rgba = [int(255*i) for i in array[:3]] + [array[3]]
    return 'rgba({:d},{:d},{:d},{:f})'.format(*int_rgba)


def combine_graphs(graphs, dynamic_load=True, **kwargs):
    # Merge graphs
    graph = None
    graphs_len = []
    for i, graph_new in enumerate(graphs):
        if dynamic_load:
            try:
                graph_new =  load_graph_by_id(graph_new, **kwargs)
            except:
                continue
        graph_new = graph_new.rename(columns={'coef': f'coef_{i}'})
        graphs_len.append(len(graph_new))
        if graph is None:
            graph = graph_new
        else:
            graph = graph.merge(
                graph_new,
                how='outer',
                on=['TF', 'TG'],
            )

    # Raise exception if no graphs
    if graph is None:
        raise LookupError('No graphs found.')

    # Replace nan with 0
    graph = graph.fillna(0)

    # Average new coefs
    graph['coef'] = graph.iloc[:, 2:].mean(axis=1)

    # Match size to originals
    # graphs_avg_len = np.mean(graphs_len)
    # graph = graph.sort_values('coef')
    # graph = graph.iloc[:int(graphs_avg_len)]

    # Cull edges for high coef
    # graph = graph[graph['coef'] > graph['coef'].quantile(.8)]

    return graph[['TF', 'TG', 'coef']]


def simulate_diffusion(g, genes, spread=.1, eps=1e-6, color_power=1):
    # Perform diffusion recursively
    g.vp.diffusion_value = g.new_vertex_property('double')
    vertices = []
    for v in g.vertices():
        if g.vp.ids[v] in genes:
            g.vp.diffusion_value[v] = 1
            # Both out and in neighbors
            vertices.append(g.get_out_neighbors(v))
            vertices.append(g.get_in_neighbors(v))
    if len(vertices) > 0:
        vertices = np.concatenate(vertices)
        g = simulate_diffusion_helper(g, vertices, 1, spread, eps)

    # Convert to color
    g.vp.diffusion_color = g.new_vertex_property('string')
    for v in g.vertices():
        # Red to white as it spreads
        val = 1-g.vp.diffusion_value[v]**color_power
        g.vp.diffusion_color[v] = rgba_to_hex((1, val, val, 1))

    return g


def simulate_diffusion_helper(g, vertices, value, spread, eps):
    if value < eps or len(vertices) == 0:
        return g
    new_value = spread * value
    new_vertices = []
    for v in vertices:
        if new_value > g.vp.diffusion_value[v]:
            g.vp.diffusion_value[v] = new_value
            new_vertices.append(g.get_out_neighbors(v))
            new_vertices.append(g.get_in_neighbors(v))
    if len(new_vertices) == 0: return g
    new_vertices = np.concatenate(new_vertices)

    return simulate_diffusion_helper(g, new_vertices, new_value, spread, eps)


def detect_synthetic_vertices_list(graph, hub_present=False):
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
    return synthetic_vertices


def detect_synthetic_vertices_graph(g):
    return [g.vp.ids[v] for v in g.vertices() if g.vp.text_synthetic[v]]


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
                indexer=lambda e: get_edge_string(gc, e),
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
            indexer=lambda e: get_edge_string(gc, e),
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
                (g.vp.shape, gc.vp.shape),
                (g.vp.ids, gc.vp.ids),
                (g.vp.text, gc.vp.text),
                (g.vp.text_synthetic, gc.vp.text_synthetic),
                (g.ep.coef, gc.ep.coef),
            ],
            include=True)
        g.vp.color, g.vp.shape, g.vp.ids, g.vp.text, g.vp.text_synthetic, g.ep.coef = props

    # Label self loops
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop)

    # Remove duplicate edges
    g = _remove_duplicate_edges(g)

    # Add processed attributes
    g.vp.self_loop_value = g.new_vertex_property('double')
    g.ep.color = g.new_edge_property('vector<double>')
    g.ep.coefs = g.new_edge_property('vector<double>')
    for e in g.edges():
        # Get coefs
        coefs = g_coefs[get_edge_string(g, e)]
        g.ep.coefs[e] = coefs

        # Get processed attributes
        in_graph = [c != 0 for c in coefs]
        present_coef = np.mean([c for c in coefs if c != 0]) if exclude_zeroes_from_mean else np.mean(coefs)
        color = [0 for _ in range(4)]

        # Set color
        color[:3] = _determine_color(g, e)

        # Set alpha
        color[3] = get_alpha(present_coef)

        # Write
        g.ep.color[e] = color
        if g.ep.self_loop:
            g.vp.self_loop_value[e.source()] = present_coef

    return g


def _determine_color(g, e, method='presence'):
    coefs = g.ep.coefs[e]

    # Get color index
    if method == 'presence':
        # Default method assigns binary string to presence
        # i.e. using two: '11' for both, '10' for second, etc.
        cindex = ''.join([str(int(b)) for b in [c != 0 for c in coefs]])[::-1]
    elif method == 'max':
        # Only assign '1' to max entry
        cindex = '0' * len(coefs)
        cindex = list(cindex)  # Dumb python string indexing workaround
        cindex[-np.argmax(coefs)] = '1'
        cindex = ''.join(cindex)
    else:
        raise AttributeError(f'Method {method} not found.')

    # Convert to color
    cindex = int(cindex, 2) # Binary to int
    if cindex == 2**len(coefs) - 1: return [0 for _ in range(3)]
    if len(coefs) == 2:
        # Custom colors
        if cindex == 1:
            hue = 2/3.
        if cindex == 2:
            hue = 0.
    else:
        hue = cindex * (1. / (2**len(coefs) - 2))
    color = colorsys.hsv_to_rgb(hue, 1, .75)

    return color


def get_edge_string(g, e):
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
    return gt.fruchterman_reingold_layout(g, weight=g.ep.coef, grid=False, scale=scale)


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


def get_alpha(coef):
    x = np.log10(1+coef)  # Log scaling
    # x = x**(1/3)  # Power scaling
    alpha = x / (1+x)
    alpha = .05 + .95 * alpha  # Add floor
    return alpha


def remove_text_by_centrality(g, percentile=95, eps=1e-10):
    # Calculate betweenness
    vertex_betweenness, _ = gt.betweenness(g)

    # Get no synthetic view
    g_nosynthetic = gt.GraphView(
        g,
        vfilt=[g.vp.text_synthetic[v] == '' for v in g.vertices()],
    )
    threshold = np.percentile([vertex_betweenness[v] for v in g_nosynthetic.vertices()], percentile)
    threshold = max(eps, threshold)  # Use eps as min threshold

    # Remove text
    for v in g_nosynthetic.vertices():
        if vertex_betweenness[v] < threshold:
            g_nosynthetic.vp.text[v] = ''

    return g


def get_intersection(g):
    return gt.GraphView(
        g,
        efilt=lambda e: sum([val > 0 for val in g.ep.coefs[e]]) == len(g.ep.coefs[e]),
    )


def cull_isolated_leaves(g):
    # Remove nodes which aren't connected to a synthetic node
    return gt.GraphView(
        g,
        vfilt=lambda v: (g.vp.text_synthetic[v] != '') or (len([n for n in v.all_neighbors() if g.vp.text_synthetic[n]]) > 0),
    )


def _remove_duplicate_edges(g):
    return gt.GraphView(
        g,
        efilt=lambda e: not _is_duplicate_edge(g, e),
    )


def _is_duplicate_edge(g, e):
    # Naive duplicate edge detection
    source = e.source()
    target = e.target()

    # Detect matches
    matches = []
    priority = False
    for f in g.edges():
        if e == f:
            # Return false if duplicate but first instance
            return False
        if source == f.source() and target == f.target():
            return True
    return False


def compute_differences(g):
    # Compute range for each edge coef
    g.ep.coef_diff = g.new_edge_property('double')
    for e in g.edges():
        coefs = [val for val in g.ep.coefs[e] if val]  # Only use nonzero values
        g.ep.coef_diff[e] = max(coefs) - min(coefs)

    return g


def color_by_significance(g):
    # Get diff if not present
    try:
        g.ep.coef_diff
    except:
        compute_differences(g)

    # Set color
    for e in g.edges():
        color = [0 for _ in range(4)]
        color[:3] = _determine_color(g, e, method='max')  # Color
        color[3] = get_alpha(g.ep.coef_diff[e])  # Opacity based on diff
        g.ep.color[e] = color

    return g


def get_default_scale(g):
    return .01 * (1 + np.log(g.num_vertices()))
