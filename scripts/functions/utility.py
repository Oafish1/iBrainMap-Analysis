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


def get_alpha(coef):
    x = np.log10(1+coef)  # Log scaling
    # x = x**(1/3)  # Power scaling
    alpha = x / (1+x)
    alpha = .05 + .95 * alpha  # Add floor
    return alpha
