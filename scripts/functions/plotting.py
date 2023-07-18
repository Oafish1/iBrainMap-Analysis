import itertools

import graph_tool.all as gt
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
import sympy

from .utility import *


### Plotting functions
def plot_tf_outgoing(graph_summary, ax=None):
    # Phenotype and Sub-Phenotype vs TF Outgoing
    sns.boxplot(
        graph_summary,
        x=graph_summary.columns[1],
        y='TF Outgoing',
        hue=graph_summary.columns[2],
        ax=ax,
    )
    if ax:
        ax.set_title('TGs per Regulon')
    else:
        plt.title('TGs per Regulon')


def plot_tf_closeness(graph_summary, ax=None):
    # Phenotype and Sub-Phenotype vs TF Outgoing
    sns.boxplot(
        graph_summary,
        x=graph_summary.columns[1],
        y='TF Closeness',
        hue=graph_summary.columns[2],
        showfliers=False,
        ax=ax,
    )


def plot_nps(meta, sample_ids, ax=None):
    data = meta
    data.index = data['SubID']

    # Filter to NPS
    nps_cols = [col for col in meta.columns if col.startswith('nps_')]
    data = meta.loc[[id in sample_ids for id in meta['SubID']], nps_cols]

    # Change values
    data = data.replace('True', True)
    data = data.replace('False', False)
    data = data.fillna(False)
    data = data.transpose()

    # Heatmap
    sns.heatmap(data=data, cbar=False, ax=ax)

    # Format
    ax.set_xticklabels(ax.get_xticklabels(), rotation=-60)


def plot_label(s, ax=None):
    ax.text(
        0, 1, s,
        horizontalalignment='right',
        verticalalignment='bottom',
        fontsize=30,
        weight='bold',
        transform=ax.transAxes)


def plot_sankey(meta, flow, order=None, use_nan=False):
    # cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Grab sankey vars
    indices = {}
    label = []
    x = []
    y = []
    label_color = []
    source = []
    target = []
    value = []
    edge_color = []

    for source_col, target_col in flow:
        # Add to labels and assign indices
        for col in (source_col, target_col):
            if col not in indices:
                indices[col] = {}
                for val in np.unique(meta[col].astype(str)):
                    if not use_nan and val == 'nan': continue
                    idx = len(label)
                    indices[col][val] = idx
                    label.append(val)
                    if order:
                        x.append(order[col])
                    y.append(len(indices[col]) - 1)
                    # label_color.append(cycle[idx])

    # Scale positions
    x = np.array(x)
    x = (x - x.min()) / (x.max() - x.min())
    x = .9 * x + .05

    y = np.array(y)
    y_new = []
    for xi in np.unique(x):
        y_sub = y[x==xi]
        y_sub = (y_sub - y_sub.min()) / (y_sub.max() - y_sub.min())
        y_sub = .9 * y_sub + .05
        y_new.append(y_sub)  # No idea why reassignment doesn't work
    y = np.concatenate(y_new)
    # Generate, shuffle, and fade
    n = sympy.nextprime(len(label))
    label_color = plt.cm.rainbow(np.linspace(0, 1, n))
    m = max(n // 3, 1)
    label_color = [label_color[i*m%len(label_color)] for i in range(len(label_color))]
    for i in range(len(label_color)):
        for j in range(3):
            label_color[i][j] = (label_color[i][j] + 1) / 2

    for source_col, target_col in flow:
        # Add edges
        for source_val, target_val in itertools.product(indices[source_col].keys(), indices[target_col].keys()):
            source.append(indices[source_col][source_val])
            target.append(indices[target_col][target_val])
            value.append(sum(
                (meta[source_col].astype(str) == source_val)
                * (meta[target_col].astype(str) == target_val)
            ))
            edge_color.append(label_color[indices[target_col][target_val]].copy())

    # Make edges translucent
    for i in range(len(edge_color)):
        edge_color[i][3] = .2

    # Convert to rgba
    label_color = [rgba_array_to_rgba_string(c) for c in label_color]
    edge_color = [rgba_array_to_rgba_string(c) for c in edge_color]

    # Plot
    fig = go.Figure(data=[go.Sankey(
        arrangement='snap',
        node=dict(
            label=label,
            # x=x,
            # y=y,
            color=label_color,
            pad=3),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=edge_color),
    )])
    # fig.update_layout(title_text='Data Overview', font_size=10)
    fig.show()


### Graph visualizations
def visualize_graph(g, scale=.06, ax=None):
    if ax is None: ax = plt.gca()
    # TODO: Prioritize synthetic nodes on top
    np.random.seed(42)
    gt.seed_rng(42)
    gt.graph_draw(
        g,
        pos=gt.sfdp_layout(g),  # sfdp_layout(g), gt.arf_layout(g, max_iter=1000), radial_tree_layout(g, root), random_layout(g)
        ink_scale=scale,
        vertex_fill_color=g.vp.color,
        vertex_size=gt.prop_to_size(g.vp.self_loop_value, 10, 30, power=1.5),
        vertex_text=g.vp.text,
        vertex_font_size=8,
        vertex_text_position=-2,  # No automatic node scaling
        edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
        edge_end_marker='arrow',
        edge_marker_size=gt.prop_to_size(g.ep.coef, 2, 7, power=1.5),
        mplfig=ax,
    )

    # Custom Legend
    palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
    legend_elements = [
        Line2D([0], [0], marker='o', color='gray', label='Hub', markerfacecolor=palette[0], markersize=15),
        Line2D([0], [0], marker='o', color='gray', label='Cell Type', markerfacecolor=palette[1], markersize=15),
        Line2D([0], [0], marker='o', color='gray', label='TF', markerfacecolor=palette[2], markersize=15),
        Line2D([0], [0], marker='o', color='gray', label='TG', markerfacecolor=palette[3], markersize=15),
        Line2D([0], [0], marker='o', color='gray', label='TF+TG', markerfacecolor=palette[4], markersize=15),
    ]
    ax.legend(handles=legend_elements, loc='best')


def visualize_graph_diffusion(g, scale=.12, ax=None):
    # TODO: Prioritize synthetic nodes on top
    np.random.seed(42)
    gt.seed_rng(42)
    gt.graph_draw(
        g,
        pos=gt.sfdp_layout(g),  # sfdp_layout(g), gt.arf_layout(g, max_iter=1000), radial_tree_layout(g, root), random_layout(g)
        ink_scale=scale,
        vertex_fill_color=g.vp.diffusion_color,
        vertex_size=gt.prop_to_size(g.vp.self_loop_value, 10, 30, power=1.5),
        vertex_text=g.vp.text,
        vertex_font_size=8,
        vertex_text_position=-2,  # No automatic node scaling
        edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
        edge_end_marker='arrow',
        edge_marker_size=gt.prop_to_size(g.ep.coef, 2, 7, power=1.5),
        mplfig=ax,
    )

    # Custom Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='gray', label='Disease-Related', markerfacecolor='#FF0000', markersize=15),
        Line2D([0], [0], marker='o', color='gray', label='Unrelated', markerfacecolor='#FFFFFF', markersize=15),
    ]
    ax.legend(handles=legend_elements, loc='best')


def visualize_graph_state(g, scale=.01, highlight=False, ax=None):
    np.random.seed(42)
    gt.seed_rng(42)
    state = gt.minimize_nested_blockmodel_dl(g)

    if highlight:
        # Gray except hub cluster
        vp_clusters = state.get_clabel(0)
        # Get hub cluster
        v_hub = [v for v in g.vertices() if g.vp.ids[v] == 'hub'][0]
        v_cluster = [v for v in g.vertices() if vp_clusters[v] == vp_clusters[v_hub]]
        # Color
        hub_cluster_color_v = g.new_vertex_property('string')
        # hub_cluster_color_e = g.new_edge_property('string')
        # TODO: Better way to default
        for v in g.vertices(): hub_cluster_color_v[v] = '#777777'
        # for e in g.edges():  hub_cluster_color_e[e] = '#777777'
        for v in g.vertices():
            if v in v_cluster:
                hub_cluster_color_v[v] = '#FF7777'
                # for e in np.concatenate([g.get_out_edges(v), g.get_in_edges(v)]):
                #     hub_cluster_color_e[e] = '#FF7777'

        state.draw(
            ink_scale=scale,
            vertex_fill_color=hub_cluster_color_v,
            # edge_color = hub_cluster_color_e,
            edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
            # vertex_text=g.vp.text_synthetic,
            vertex_text_position='centered',
            vertex_font_size=6,
            mplfig=ax,
        )

    else:
        state.draw(
            ink_scale=scale,
            edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
            mplfig=ax,
        )
