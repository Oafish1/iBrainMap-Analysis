import itertools

import graph_tool.all as gt
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
import sympy

from .computation import *
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


def plot_statistic(graph_summary, col='TF Closeness', ax=None):
    # Phenotype and Sub-Phenotype vs TF Outgoing
    sns.boxplot(
        graph_summary,
        x=graph_summary.columns[1],
        y=col,
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
    sns.heatmap(data=data, square=True, xticklabels=True, yticklabels=True, cbar=False, ax=ax)

    # Format
    ax.set_xticklabels(ax.get_xticklabels(), rotation=-90)


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

    if order:
        # Scale positions
        x = np.array(x)
        x = (x - x.min()) / (x.max() - x.min())
        x = .9 * x + .05

        y = np.array(y)
        y_new = []
        for xi in np.unique(x):
            y_sub = y[x==xi]
            y_sub = (y_sub - y_sub.min()) / (y_sub.max() - y_sub.min())
            y_sub = .8 * y_sub + .1  # Lower is higher
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
            pad=15,
            thickness=15,
            label=label,
            x=x if order else None,
            y=y if order else None,
            color=label_color),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=edge_color),
    )])
    fig.update_layout(title_text='Data Overview', font_size=10)
    return fig


### Graph visualizations
def visualize_graph(g, pos=None, scale=None, ax=None, legend=False):
    if not scale: scale = get_default_scale(g)
    if not ax: ax = plt.gca()
    if not pos: pos = get_graph_pos(g)
    # TODO: Prioritize synthetic nodes on top
    np.random.seed(42)
    gt.seed_rng(42)
    gt.graph_draw(
        g,
        pos=pos,
        ink_scale=scale,
        vertex_fill_color=g.vp.color,
        vertex_shape=g.vp.shape,
        vertex_size=gt.prop_to_size(g.vp.self_loop_value, 20, 40, power=1.5),
        vertex_text=g.vp.text,
        vertex_font_size=8,
        vertex_text_position=-2,  # No automatic node scaling
        # edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
        edge_color=g.ep.color,
        edge_end_marker='arrow',
        edge_marker_size=10,  # gt.prop_to_size(g.ep.coef, 2, 7, power=1.5),
        mplfig=ax,
    )

    if legend: plot_legend


def plot_legend(hub=False, ax=None):
    if not ax: ax = plt.gca()

    # Custom Legend
    palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
    legend_elements = [
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='p', markerfacecolor=palette[1], label='Cell Type'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=palette[4], label='TF+TG'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=palette[2], label='TF'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='o', markerfacecolor=palette[3], label='TG'),
    ]
    if hub:
        legend_elements = [Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='h', markerfacecolor=palette[1], label='Hub'),] + legend_elements
    ax.legend(handles=legend_elements, loc='best')


def visualize_graph_diffusion(g, pos=None, scale=None, ax=None):
    if not scale: scale = get_default_scale(g)
    if not ax: ax = plt.gca()
    if not pos: pos = get_graph_pos(g)
    # TODO: Prioritize synthetic nodes on top
    np.random.seed(42)
    gt.seed_rng(42)
    gt.graph_draw(
        g,
        pos=pos,
        ink_scale=scale,
        vertex_fill_color=g.vp.diffusion_color,
        vertex_shape=g.vp.shape,
        vertex_size=gt.prop_to_size(g.vp.self_loop_value, 20, 40, power=1.5),
        vertex_text=g.vp.text,
        vertex_font_size=4,
        vertex_text_position=-2,  # No automatic node scaling
        # edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
        edge_color=g.ep.color,
        edge_end_marker='arrow',
        edge_marker_size=10,  # gt.prop_to_size(g.ep.coef, 2, 7, power=1.5),
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


def plot_enrichment(df, ax=None):
    if not ax: ax = plt.gca()

    # Plot
    sns.scatterplot(
        data=df,
        x='Disease',
        y='Cell Type',
        size='-log10(p)',
        color='black',
        sizes=(10, 200),
        ax=ax,
    )
    # Zoom X
    margin = .5
    min_xlim, max_xlim = ax.get_xlim()
    min_xlim -= margin; max_xlim += margin
    # Zoom Y
    max_ylim, min_ylim = ax.get_ylim()
    min_ylim -= margin; max_ylim += margin
    ax.set(xlim=(min_xlim, max_xlim), ylim=(max_ylim, min_ylim))
    # Formatting
    plt.grid()
    plt.xticks(rotation=90)
    # Legend
    plt.legend(bbox_to_anchor=(1.02, .7), loc='upper left', borderaxespad=0, frameon=False)  # Legend to middle-right outside


def plot_individual_edge_comparison(g, sample_ids, ax=None):
    if not ax: ax = plt.gca()

    # Assemble dataframe
    df = {k: [] for k in ['id'] + sample_ids}
    for e in g.edges():
        coefs = g.ep.coefs[e]
        df['id'].append(get_edge_string(g, e))
        for i, sample_id in enumerate(sample_ids):
            df[sample_id].append(coefs[i])
    df = pd.DataFrame(df)

    # Plot
    sns.scatterplot(
        data=df,
        x=sample_ids[0],
        y=sample_ids[1],
        color='black',
        ax=ax,
    )

    # Plot y=x
    lims = [
            max(ax.get_xlim()[0], ax.get_ylim()[0]),
            min(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, '-', color='black', alpha=0.3)

    # Formatting
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    return df
