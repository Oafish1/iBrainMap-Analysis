
import colorsys
import itertools
import os

from adjustText import adjust_text
from brokenaxes import brokenaxes
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
    "Plotting with additional parameters, unused"
    if not scale: scale = get_default_scale(g)
    if not ax: ax = plt.gca()
    if not pos: pos = get_graph_pos(g)
    # TODO: Prioritize synthetic nodes on top
    np.random.seed(42)
    gt.seed_rng(42)
    visualize_graph_base(
        g,
        pos=pos,
        ink_scale=scale,
        vertex_size=gt.prop_to_size(g.vp.self_loop_value, 20, 40, power=1.5),
        vertex_font_size=8,
        vertex_text_position=-2,  # No automatic node scaling
        # edge_pen_width=gt.prop_to_size(g.ep.coef, .1, 1, power=1.5),
        edge_marker_size=10,  # gt.prop_to_size(g.ep.coef, 2, 7, power=1.5),
        mplfig=ax,
    )

    if legend: plot_legend()


def visualize_graph_base(g, **kwargs):  # , min_size=None
    "Basic graph visualization"
    # Scale size
    # if min_size is None: min_size = 1*g.num_vertices()**(-5/6)
    min_size = .7*g.num_vertices()**(-4/6)
    max_size = 5*min_size
    size = g.new_vertex_property('double')
    for v in g.vertices():
        size[v] = max_size * g.vp.size[v] + min_size

    # Only show visible nodes
    # This method is strange, should revise
    outline_color = g.new_vertex_property('vector<double>')
    for v in g.vertices():
        try:
            if not g.vp.hide[v]: outline_color[v] = [0, 0, 0, .2]
        except: outline_color[v] = [0, 0, 0, .2]


    # Draw
    gt.graph_draw(
        g,
        vertex_fill_color=g.vp.color,
        vertex_color=outline_color,
        vertex_shape=g.vp.shape,
        vertex_size=size,  # gt.prop_to_size(g.vp.size, min_size, max_size),  # vertex_size uses absolute units
        vertex_text=g.vp.text,
        vertex_text_position=-2,
        vertex_text_color='black',
        edge_color=g.ep.color,
        edge_end_marker='none',
        # ink_scale=1,
        fit_view=1.1,
        **kwargs,
    )

    # return min_size


def plot_legend(horizontal=True, hub=False, loc='best', ax=None, **kwargs):
    if not ax: ax = plt.gca()

    # Custom Legend
    # NOTE: Unfortunately, shapes need to be changed manually
    legend_elements = [
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='h', markerfacecolor=get_node_appearance('celltype')[0], label='Cell Type'),
        # Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=get_node_appearance('tftg')[0], label='TF+TG'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=get_node_appearance('tf')[0], label='Transcription Factor'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='o', markerfacecolor=get_node_appearance('tg')[0], label='Target Gene'),
    ]
    if hub:
        legend_elements = [Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='8', markerfacecolor=get_node_appearance('hub')[0], label='Hub'),] + legend_elements
    ax.legend(handles=legend_elements, loc=loc, ncol=1 if not horizontal else len(legend_elements), **kwargs)


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
    "Macro for `plot_circle_heatmap` for enrichment results"
    plot_circle_heatmap(
        df,
        index_name='Cell Type',
        column_name='Disease',
        value_name='-log(p)',
        color='Black',
        transform=False,
        ax=ax)
    ax.set_xlabel(None)
    ax.set_ylabel(None)


def plot_individual_edge_comparison(g, sample_ids, suffix='Prioritization Weights', highlight_outliers=True, color_map=None, filter=lambda x: filter_remove_ct_ct(filter_remove_tf_tg(x)), ax=None):
    """
    Take concatenated graph `g` and plot a comparison between the original weights

    By default, only plots CT-TF edges
    """
    if not ax: ax = plt.gca()

    # Assemble dataframe
    df = {k: [] for k in ['id'] + sample_ids}
    for e in g.edges():
        coefs = g.ep.coefs[e]
        df['id'].append(get_edge_string(g, e))
        for i, sample_id in enumerate(sample_ids):
            df[sample_id].append(coefs[i])
    df = pd.DataFrame(df)

    # Filter
    if filter is not None: df = filter(df).reset_index(drop=True)

    # Highlight outliers
    if highlight_outliers:
        # Get outliers
        difference = (df[sample_ids[1]] - df[sample_ids[0]]).to_numpy()
        outlier_idx, outlier_mask = get_outlier_idx(difference, return_mask=True)
        # Format based on which is higher
        # 2 if sid1 is higher, 1 otherwise
        sid_mask = 1. * (difference > 0) + 1
        # Record
        df['outlier'] = sid_mask * outlier_mask
        # Annotate
        annotations = []
        for idx in outlier_idx:
            # `ax.annotate` doesn't work with `adjust_text`
            # upper = difference[idx] > 0  # Is the point above or below y=x?
            # ann = ax.annotate(
            #     df.loc[idx, 'id'].item(),
            #     (df.loc[idx, sample_ids[0]].item(), df.loc[idx, sample_ids[1]].item()),
            #     # Away from y=x
            #     # horizontalalignment='right' if upper else 'left',
            #     # verticalalignment='bottom' if upper else 'top',
            #     horizontalalignment='center',
            #     verticalalignment='center',
            #     # Formatting
            #     fontsize=10,
            # )
            ann = ax.text(
                df.loc[idx, sample_ids[0]].item(), df.loc[idx, sample_ids[1]].item(),
                df.loc[idx, 'id'].item(),
                # Away from y=x
                # horizontalalignment='right' if upper else 'left',
                # verticalalignment='bottom' if upper else 'top',
                ha='center',
                va='center',
                # Formatting
                fontsize=10,
            )
            annotations.append(ann)

    # Plot
    sns.scatterplot(
        data=df,
        x=sample_ids[0],
        y=sample_ids[1],
        hue='outlier' if highlight_outliers else None,
        palette={
            0: 'black',
            1: color_map[sample_ids[0]] if color_map is not None else 'red',
            2: color_map[sample_ids[1]] if color_map is not None else 'red'},
        color='black',
        # linewidth=0,
        ax=ax,
    )
    xlabel = sample_ids[0]
    if suffix: xlabel += ' ' + suffix
    ax.set_xlabel(xlabel)
    ylabel = sample_ids[1]
    if suffix: ylabel += ' ' + suffix
    ax.set_ylabel(ylabel)
    plot_remove_legend(ax=ax)

    # Formatting
    # NOTE: Do before y=x for tighter boundaries
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_aspect('equal', 'box')

    # Plot y=x
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    lims = [
        max(ax.get_xlim()[0], ax.get_ylim()[0]),
        min(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, '-', color='black', alpha=0.3)

    # Adjust text positions to avoid overlap
    # NOTE: Must be done after axis limit/scale changes
    if highlight_outliers and len(annotations) > 0:
        adjust_text(
            annotations,
            df[sample_ids[0]].to_numpy(),
            df[sample_ids[1]].to_numpy(),
            arrowprops=dict(arrowstyle='-', color='black', lw=.01),
            ax=ax,
        )

    return df


def get_mosaic(mosaic, figsize=None, **kwargs):
    """
    Produce a mosaic using `mosaic` as a layout

    figsize: (length, height)
    """
    if figsize is None: figsize = (6, 6)
    fig = plt.figure(figsize=figsize, **kwargs)
    axs = fig.subplot_mosaic(mosaic)

    return fig, axs


def plot_graph_comparison(
        graphs=None,
        *,
        concatenated_graph=None,
        concatenated_pos=None,
        axs,
        subject_ids,
        filter_text=False,
        show_null_nodes=True,
        include_titles=True,
        **kwargs):
    # Aggregate and calculate node positions
    # NOTE: Supplied concatenated graph must not be pruned
    if concatenated_graph is None:
        concatenated_graph = concatenate_graphs(*graphs, recalculate=False, threshold=False)
        if filter_text:
            concatenated_graph = remove_text_by_centrality(concatenated_graph.copy())
    if concatenated_pos is None:
        concatenated_pos = get_graph_pos(concatenated_graph)

    # Plot
    for i, (g, sid) in enumerate(zip(graphs, subject_ids)):
        ax = axs[i]
        if i:
            axs[0].sharex(ax)
            axs[0].sharey(ax)
        if show_null_nodes:
            inverse_graph = (
                make_vertices_white(
                make_snps_invisible(
                get_inverse_graph(
                remove_edges(
                remove_text(
                    concatenated_graph.copy()
                )), g)))
            )
            plot_graph = concatenate_graphs(
                inverse_graph,
                transfer_text_labels(concatenated_graph, g),
                recolor=False,
                threshold=False,
                recalculate=False)
        else: plot_graph = transfer_text_labels(concatenated_graph, g)
        visualize_graph_base(
            plot_graph,
            pos=convert_vertex_map(concatenated_graph, plot_graph, concatenated_pos),
            # vertex_font_size=.5*20/g.num_vertices()**(1/2),
            mplfig=ax,
            **kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        if include_titles: ax.set_title(sid)
        ax.axis('off')


def plot_edge_summary(graphs, *, df=None, ax, subject_ids=None, min_common_edges=1, num_x_labels=15):
    if df is None:
        assert subject_ids is not None
        df, _ = compute_edge_summary(graphs=graphs, subject_ids=subject_ids, min_common_edges=min_common_edges)

    # Melt by subject ID
    df_long = pd.melt(df, id_vars=['Edge', 'Variance', 'Mean'], var_name='Subject ID', value_name='Weight')
    df_long = df_long.sort_values(['Variance'], ascending=[True])

    # Don't plot zero values
    df_long = df_long.loc[df_long['Weight']!=0]

    # Plot
    lp = sns.lineplot(
        data=df_long,
        x='Edge', y='Weight', hue='Subject ID',
        alpha=.5,
        errorbar=None, estimator=None, n_boot=0,
        ax=ax)
    plt.xlabel(None)
    plt.xticks(rotation=90)

    # Only show `num_x_labels` x labels
    limit_labels(lp, n=num_x_labels)


def plot_aggregate_edge_summary(contrast_subject_ids=None, *, contrast=None, column=None, ax, min_common_edges=1, num_x_labels=15, **kwargs):
    # Setup
    if contrast is not None: contrast_concatenated_graphs, contrast_concatenated_subject_ids = contrast
    else:
        assert column is not None
        contrast_concatenated_graphs, contrast_concatenated_subject_ids = compute_aggregate_edge_summary(contrast_subject_ids, column=column, **kwargs)

    # Aggregate data
    aggregate_df = pd.DataFrame(columns=['Edge', 'Subgroup', 'Variance', 'Mean'])
    for key, concatenated_graph in contrast_concatenated_graphs.items():
        df, _ = compute_edge_summary(concatenated_graph=concatenated_graph, subject_ids=contrast_concatenated_subject_ids[key])
        df['Subgroup'] = key
        aggregate_df = pd.concat([aggregate_df, df[list(aggregate_df.columns)]])

    # Remove edges with few readings
    unique, counts = np.unique(aggregate_df['Edge'], return_counts=True)
    aggregate_df = aggregate_df.loc[aggregate_df['Edge'].isin(unique[counts >= min_common_edges])]

    # Get mean difference
    # TODO: Plot mean difference

    # Sort
    aggregate_df = aggregate_df.sort_values(['Variance'], ascending=[True])

    # Plot
    lp = sns.lineplot(  # Takes a while
        data=aggregate_df,
        x='Edge', y='Mean', hue='Subgroup',
        alpha=.5,
        errorbar=None, estimator=None, n_boot=0,
        ax=ax)
    plt.xlabel(None)
    plt.xticks(rotation=90)

    # Only show `num_x_labels` x labels
    len_loop = len(lp.get_xticklabels())
    for i, label in enumerate(lp.get_xticklabels()):
        if i % int(len_loop/num_x_labels) != 0: label.set_visible(False)

    return aggregate_df


def get_graph_pos(g, scale=None):
    "Make scaling consistent for positioning in graph_tools"
    # Compute scale
    # if not scale:
    #     scale = 20 * np.sqrt(g.num_vertices())

    # Compute layout
    print('Calculating positions...')
    # sfdp_layout(g), gt.arf_layout(g, max_iter=1000), radial_tree_layout(g, root), random_layout(g)
    # return gt.random_layout(g)
    # return gt.arf_layout(g, weight=g.ep.coef)

    # Random
    # pos = gt.random_layout(g)

    # SFDP
    pos = gt.sfdp_layout(
        g,
        # vweight=g.vp.self_loop_value,
        # eweight=g.ep.coef,
        # kappa=.5, gamma=.1, r=2,
        # r=50,
        # p=4,  # Fewer stragglers very far away
    )

    # Fruchterman Reingold
    # pos = gt.fruchterman_reingold_layout(
    #     g,
    #     # weight=g.ep.coef,
    #     # grid=False,
    #     # scale=scale,
    #     # n_iter=100,
    # )

    # Debug
    # for v in g.vertices():
    #     print(pos[v])

    pos = scale_pos_to_range(g, pos)
    pos = scale_pos_by_distance(g, pos)
    pos = scale_pos_to_range(g, pos)
    return pos


def plot_contrast_curve(
    df_subgroup: pd.DataFrame,
    *,
    index_name: str='Edge',
    subgroup_name: str='Subgroup',
    value_name: str='Variance',
    sorting_subgroup: str='Individual',
    concatenate: bool=True,
    filter_common: bool=None,
    ax=None,
    **kwargs):
    """
    Plot the variance curve for subgroups in `df_subgroup`.

    `sorting_subgroup`: Subgroup to sort by.  If 'Individual',
      Sort each subgroup on its own.
    `filter_common`: Filter to common edges between all subgroups.
      Necessary for comparability.
    """
    # TODO: Add filter_super option which filters only to sorting group
    # Defaults
    if filter_common is None: filter_common = sorting_subgroup != 'Population'
    if ax is None: _, axs = get_mosaic([list(range(1))], scale=9); ax=axs[0]

    # Merge and format dfs
    if concatenate:
        df_concat = pd.concat([df_subgroup[k] for k in df_subgroup])
    else:
        df_concat = df_subgroup
    df_concat = df_concat.set_index(index_name)

    # Filter to only common edges if sorting subgroup is not population, for comparability
    if filter_common:
        unique, count = np.unique(df_concat.index, return_counts=True)
        edges_to_keep = unique[count >= len(np.unique(df_concat[subgroup_name]))]
        df_concat = df_concat.loc[[s in edges_to_keep for s in df_concat.index]]

    # Sort
    if sorting_subgroup == 'Individual':
        # Have each line individually sorted
        df_concat = df_concat.sort_values([value_name], ascending=[True])
        df_concat['Sort'] = 0
        for subgroup in np.unique(df_concat[subgroup_name]):
            length = len(df_concat.loc[df_concat[subgroup_name] == subgroup, 'Sort'])
            df_concat.loc[df_concat[subgroup_name] == subgroup, 'Sort'] = (
                np.array(list(range(length))) / (length - 1))
        x = 'Sort'
        xlabel = 'Percentile'
    elif sorting_subgroup == 'Mean':
        # Sort by mean of all except population
        df_concat_edge_mean = (
            df_concat
                .loc[df_concat[subgroup_name] != 'Population', [value_name]]
                .groupby([index_name])
                .mean()
        )
        df_concat = df_concat.join(df_concat_edge_mean, how='right', rsuffix='_sort')
        df_concat = df_concat.sort_values([f'{value_name}_sort'], ascending=[True])
        x = index_name
        xlabel = None
    else:
        # Sort by a single column
        # Prepare for join
        to_join = df_concat.loc[df_concat[subgroup_name] == sorting_subgroup]
        # to_join.index = to_join[index_name]
        # Join and sort
        df_concat = df_concat.join(to_join, how='right', rsuffix='_sort')  # Filter to and sort by subgroup
        df_concat = df_concat.sort_values([f'{value_name}_sort'], ascending=[True])
        x = index_name
        xlabel = None

    # Plot
    plt.sca(ax)
    lp = sns.lineplot(
        data=df_concat,
        x=x, y=value_name, hue=subgroup_name,
        hue_order=np.unique(df_concat[subgroup_name]),
        alpha=.65,
        errorbar=None, estimator=None, n_boot=0,
        ax=ax,
        **kwargs)
    plt.yscale('log')
    plt.xlabel(xlabel)
    plt.xticks(rotation=90)

    # Only show `num_x_labels` x labels
    num_x_labels = 20
    len_loop = len(lp.get_xticklabels())
    if len_loop > num_x_labels:
        for i, label in enumerate(lp.get_xticklabels()):
            if i % int(len_loop/num_x_labels) != 0: label.set_visible(False)


def plot_subgroup_heatmap(df_subgroup, *, ax=None):
    # Regular heatmap
    # sns.heatmap(data=join_df_subgroup(df_subgroup), ax=ax)

    # Circle heatmap
    plot_circle_heatmap(join_df_subgroup(df_subgroup), ax=ax)
    ax.set_xlabel(None)
    ax.set_ylabel(None)


def plot_BRAAK_comparison(contrast, *, meta, column, target='BRAAK_AD', df=None, legend=True, ax=None, **kwargs):
    # TODO: Rename?  It can do more than BRAAK
    # Compute
    if df is None:
        df, edges_include = compute_BRAAK_comparison(contrast, meta=meta, column=column, target=target, **kwargs)

    # Plot
    sns.violinplot(data=df, hue=target, y='Attention', x='Edge', ax=ax)
    # plt.yscale('log')  # Could this misfire?
    # sns.despine(offset=10, ax=ax)  # trim=True
    if legend: plot_outside_legend(title=target, ax=ax)
    else: plot_remove_legend(ax=ax)

    return df, edges_include


def plot_prediction_confusion(
        contrast,
        *,
        meta,
        column,
        target='BRAAK_AD',
        prioritized_edges,
        plot_type='Barplot',
        row_normalize: bool = None,
        ax=None,
        **kwargs):
    # Parameters
    if ax is None: ax = plt.gca()
    if row_normalize is None:
        row_normalize = plot_type.upper() == 'HEATMAP'

    # Compute
    df, acc = compute_prediction_confusion(contrast, meta=meta, column=column, target=target, prioritized_edges=prioritized_edges, **kwargs)

    # Get num samples
    n = df.to_numpy().sum()

    # Row scale df
    if row_normalize:
        df = ( df.T / df.sum(axis=1) ).T  # Kind of hacky, works because columns are default divide

    # Plot
    if plot_type.upper() == 'HEATMAP':
        sns.heatmap(data=df, vmin=0, cmap='crest', cbar=not row_normalize, ax=ax)
        ax.set_title(f'n={n}, acc={acc:.3f}')
        ax.set_xlabel(f'{target} (Predicted)')
        ax.set_ylabel(f'{target} (True)')

    elif plot_type.upper() == 'BARPLOT':
        # df_bar = df.reset_index(names='True').melt(
        #     id_vars='True',
        #     value_vars=df.columns,
        #     value_name='Count',
        #     var_name='Predicted')
        df.plot(kind='bar', stacked=True, ax=ax)
        ax.set_title(f'n={n}, acc={acc:.3f}')
        ax.set_xlabel(f'{target}')
        ax.set_ylabel(f'Count')



    return df, acc


def plot_circle_heatmap(
    df,
    *,
    index_name=None,
    column_name='Group',
    value_name='Variance',
    sign_name='Sign',
    color='Red',
    size_max=200,
    multicolor=None,
    multicolor_labels=['Positive', 'Negative'],
    transform=True,
    ax=None,
    **kwargs):
    "Plot a heatmap with circles, like some corrplots in R"
    # TODO: Set up custom colors with multicolor
    # For color, generally, Significance=Black, Variance=Red, Raw Difference=Blue
    # Setup
    if ax is None: ax = plt.gca()
    else: plt.sca(ax)
    if transform and multicolor is None: multicolor = df.to_numpy().min() < 0
    if transform and multicolor is None: multicolor = df.to_numpy().min() < 0

    # Melt and format heatmap df (like `sns.heatmap` input) for scatter
    if transform:
        # Make copies
        df = df.copy()

        # Format data
        if index_name is None: index_name = df.index.name
        df = df.reset_index(names=index_name)
        df = pd.melt(df, id_vars=[index_name], var_name=column_name, value_name=value_name)

    ## Get metadata
    # Split by positive/negative
    if multicolor:
        df[sign_name] = list(df[value_name].apply(lambda x: multicolor_labels[0] if x > 0 else multicolor_labels[1]))
        df[value_name] = df[value_name].abs()

    # Scale from zero
    size_min = size_max * df[value_name].min() / df[value_name].max()

    # Plot
    # TODO: Set min size to represent true 0, get constant scaling
    scp = sns.scatterplot(
        data=df,
        x=column_name,
        y=index_name,
        hue=sign_name if multicolor else None,
        size=value_name,
        sizes=(size_min, size_max),
        color=color,
        ax=ax,
        **kwargs,
    )
    # Formatting
    plt.grid()
    plt.xticks(rotation=90)
    h, l = scp.axes.get_legend_handles_labels()
    scp.axes.legend_.remove()
    plt.legend(h, l, ncol=2)
    # ax.set_aspect('equal', 'box')
    # Zoom X
    margin = 1
    min_xlim, max_xlim = ax.get_xlim()
    min_xlim -= margin; max_xlim += margin
    # Zoom Y
    max_ylim, min_ylim = ax.get_ylim()
    min_ylim -= margin; max_ylim += margin
    ax.set(xlim=(min_xlim, max_xlim), ylim=(max_ylim, min_ylim))
    # Legend
    plot_outside_legend(title=value_name if not multicolor else None, ax=ax)


def plot_outside_legend(title=None, frameon=False, ax=None):
    "Plot a legend to the right of plot"
    # TODO: Add title to legend
    if ax is not None: plt.sca(ax)
    plt.legend(title=title, bbox_to_anchor=(1.02, .7), loc='upper left', borderaxespad=0, frameon=frameon)  # Legend to middle-right outside


def plot_remove_legend(ax=None):
    "Remove legend on the provided plot"
    if ax is not None: plt.sca(ax)
    plt.legend([],[], frameon=False)


def plot_head_comparison(subject_id_1, subject_id_2, *, colors=None, ax=None, **kwargs):
    # Setup
    if ax is None: ax = plt.gca()

    # Compute
    individual_head_difference = compute_head_comparison((subject_id_1, subject_id_2), **kwargs)

    # Plot
    plot_circle_heatmap(
        individual_head_difference.T,
        index_name='Head',
        column_name='Edge',
        value_name='Difference',
        sign_name='Subject',
        color='Blue',  # Fallback color
        multicolor_labels=[subject_id_1, subject_id_2],
        ax=ax,
        hue_order=[subject_id_1, subject_id_2],
        palette=colors,
    )
    ax.set_xlabel(None)
    ax.set_ylabel(None)


def plot_attention_dosage_correlation(
        dosage,
        *,
        n=None,  # Helps with development for shorter runtimes
        random_state=42,
        target_edge=None,
        chromosomes=None,
        ax=None,
        **kwargs):
    # Parameters
    if ax is None: ax = plt.gca()

    # Compute correlations
    ret = compute_attention_dosage_correlation(
        dosage if n is None else dosage.sample(n=n, random_state=random_state),
        target_edge=target_edge,
        chromosomes=chromosomes,
        **kwargs)
    if target_edge is None:
        df, target_edge = ret
    else:
        df = ret

    # Plot
    plt.sca(ax)
    ax.axhline(y=-np.log(.05))  # Significance line
    sns.scatterplot(
        data=df,
        x='Genomic Coordinate (bp)',
        y='-log(FDR-Adjusted Correlation p-value)',
        hue='Chromosome',
        hue_order=['chr'+s for s in (get_chromosome_order() if chromosomes is None else chromosomes)],
        edgecolor='none',
        ax=ax)
    title = target_edge
    if chromosomes is not None: title += ' ' + ', '.join(chromosomes)
    ax.set_title(target_edge)
    if chromosomes is None: plot_remove_legend(ax=ax)

    return ret


def plot_graph_comparison_from_sids(sids, *, axs, column=None, vertex_ids=None, include_titles=False):
    """
    axs: Dictionary of axes with indices 0 and 1
    """

    # Assemble
    gs = [compute_graph(load_graph_by_id(sid, column=column)) for sid in sids]

    # Filter
    if vertex_ids is not None:
        gs = [
            filter_to_synthetic_vertices(g.copy(), vertex_ids=vertex_ids)
            for g in gs
        ]

    # Recalculate
    gs = [assign_vertex_properties(g) for g in gs]

    # Plot
    plot_graph_comparison(gs, axs=axs, subject_ids=sids, include_titles=include_titles)


def plot_edge_comparison_from_sids(sids, *, ax, column=None, palette=None, **kwargs):
    # Assemble
    graphs = [compute_graph(load_graph_by_id(sid, column=column)) for sid in sids]

    # Get graph
    g = concatenate_graphs(*graphs, threshold=False)
    g = get_intersection(g)
    # g = cull_isolated_leaves(g)

    plot_individual_edge_comparison(g, sids, color_map=palette, ax=ax, **kwargs)


def plot_module_scores_from_sids(sids, *, ax, palette, column=None):
    # Parameters
    subject_id_1, subject_id_2 = sids

    # Get graphs
    g1 = compute_graph(load_graph_by_id(subject_id_1, column=column, average=True))
    g2 = compute_graph(load_graph_by_id(subject_id_2, column=column, average=True))

    # Get module scores
    module_scores_1 = get_module_scores(g1)
    module_scores_2 = get_module_scores(g2)
    # Make blanks
    zeros_1 = module_scores_1.copy()
    zeros_1['Module Score'] = 0
    zeros_2 = module_scores_2.copy()
    zeros_2['Module Score'] = 0
    # Append for consistency
    module_scores_1 = pd.concat((module_scores_1, zeros_2)).groupby(['Cell Type', 'TF']).max().reset_index()
    module_scores_2 = pd.concat((module_scores_2, zeros_1)).groupby(['Cell Type', 'TF']).max().reset_index()
    # Concatenate subjects
    # NOTE: Only matters that they're in the order sub_1 -> sub_2
    # and all present for the `.diff()` groupby, no need to label
    # module_scores_1['Subject'] = subject_id_1
    # module_scores_2['Subject'] = subject_id_2
    module_scores = pd.concat((module_scores_1, module_scores_2))
    module_scores['Module Score'] = -module_scores.groupby(['Cell Type', 'TF'])['Module Score'].diff(periods=-1)  # Second minus first
    module_scores = module_scores.loc[~module_scores['Module Score'].isna()]

    # Filter to only high values
    # module_scores = module_scores.loc[module_scores['Module Score'].abs() > .1]

    # Plot
    def plot_module_scores(module_scores, ax=None):
        # Pivot
        df = module_scores.pivot(index='Cell Type', columns='TF', values='Module Score')
        # Roughly sort by cell type
        df = df.T
        for c in df.columns:
            df = df.sort_values(c)
        df = df.T  # .iloc[::-1]
        # Plot
        from matplotlib.colors import SymLogNorm
        pl = sns.heatmap(
            data=df,
            vmin=np.abs(df.fillna(0).to_numpy()).max(),
            vmax=-np.abs(df.fillna(0).to_numpy()).max(),
            norm=SymLogNorm(linthresh=1),
            cmap=sns.diverging_palette(
                360*colorsys.rgb_to_hls(*hex_to_rgb(palette[subject_id_1]))[0],
                360*colorsys.rgb_to_hls(*hex_to_rgb(palette[subject_id_2]))[0],
                s=80,
                l=70,
                center='dark',
                as_cmap=True),
            cbar_kws={'label': f'{subject_id_1} (Prioritization) {subject_id_2}'},
            ax=ax)
        pl.set(xlabel=f'TF (+) (n={df.shape[0]})')
        return pl
    p1 = plot_module_scores(module_scores, ax=ax)

    # Inset axis
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(
        ax,
        width='30%', height='15%',
        loc=4,
        bbox_to_anchor=(.05, .05, 1, 1), bbox_transform=ax.transAxes)
    # Take absolute module score for histogram, unevenly distributed
    # TODO: Use symlog
    module_scores_abs = module_scores.copy()
    module_scores_abs['Module Score'] = module_scores_abs['Module Score'].abs()
    sns.histplot(data=module_scores, x='Module Score', log_scale=True, kde=True, ax=axins)
    plt.ylabel(None)
    plt.xlabel(None)

    # Format
    p1.set(title=column)


def plot_edge_discovery(
        data,
        *,
        edges,
        heads,
        ax,
        column=None,
        percentage_prioritizations_ranges=[(center-center/10, center+center/10) for center in (.01, .05, .1)],
        num_labels=15,
        **kwargs):
    # Defaults
    if column is None: column = get_attention_columns()[0]

    # Get all compatible edges
    counts = compute_edge_counts(data=data, edges=edges, heads=heads)

    # Plot for chosen head
    # Filter to column
    counts_filtered = counts.loc[counts['Head']==column]

    # Sample
    # NOTE: Maybe remove in final version?  Doesn't matter too much
    np.random.seed(42)
    idx = np.random.choice(counts_filtered.shape[0], min(1_000, counts_filtered.shape[0]), replace=False)
    counts_filtered = counts_filtered.iloc[idx]

    # Sort
    counts_filtered = counts_filtered.sort_values('Count')

    # Plot
    pl = sns.lineplot(data=counts_filtered, x='Edge', y='Count', ax=ax)

    # Highlight area
    for low, high in percentage_prioritizations_ranges:
        ax.axhspan(
            low*data.shape[2],
            high*data.shape[2],
            color='red', alpha=.2, lw=0)

    # Format
    plt.sca(ax)
    plt.xticks(rotation=90)
    pl.set(title=column)
    # plt.yscale('log')
    limit_labels(pl, n=num_labels)
    # NOTE: I have utterly no idea why, but this is necessary for
    # constrained layout to work with >100 data points and not be
    # absurdly long in subplots
    # ax.set_ylabel(None)


def plot_enrichment_from_fname(fname, *, ax, num_descriptors=10):
    enrichment = pd.read_csv(fname)

    # Format
    enrichment = format_enrichment(enrichment, filter=num_descriptors)

    # Filter to certain groups
    # TODO: Fix these groups in previous section
    gene_sets = np.unique(enrichment['Gene Set'])[[0, 1, 2, 3, 4, 5, 6, 8, 10]]
    enrichment = enrichment.loc[enrichment['Gene Set'].apply(lambda s: s in gene_sets)]

    # Plot
    pl = sns.scatterplot(
        enrichment,
        x='Description', y='Gene Set',
        size='-log10(p)',
        color='black',
        ax=ax)
    # Formatting
    pl.grid()
    plt.sca(ax)
    plt.xticks(rotation=90)
    # pl.set_aspect('equal', 'box')
    pl.legend(title='-log10(p)', loc='upper left', bbox_to_anchor=(1.05, 1))
    # Zoom X
    margin = .5
    min_xlim, max_xlim = pl.get_xlim()
    min_xlim -= margin; max_xlim += margin
    pl.set(xlim=(min_xlim, max_xlim))


def plot_ct_edges(contrast, *, ax, data, edges, heads, subject_ids, column=None, num_edges=10):
    # Defaults
    if column is None: column = get_attention_columns()[0]

    # Get contrast sids
    contrast_subjects = get_contrast(contrast)

    # Assert only two groups
    assert len(contrast_subjects) == 2, 'Currently, only contrasts with two groups are supported'

    # Filter to column
    data = data[:, np.array(heads)==column]

    # Filter to ct-ct edges
    ct_edges_mask = [string_is_synthetic(s.split(get_edge_string())[0]) and string_is_synthetic(s.split(get_edge_string())[1]) for s in edges]
    data = data[ct_edges_mask]
    ct_edges = np.array(edges)[ct_edges_mask]

    # Filter to edges that have over 10 entries
    data = data[((~np.isnan(data)).sum(axis=2).flatten() > 10).astype(bool)]

    # Find means across contrast groups
    group_means = {}
    for group, sids in contrast_subjects.items():
        filtered_data = data[:, :, [sid in sids for sid in subject_ids]]
        group_means[group] = np.nanmean(filtered_data, axis=2).flatten()

    # Compute difference between means
    difference = np.abs(group_means[list(group_means.keys())[0]] - group_means[list(group_means.keys())[1]])

    # Sort and filter
    difference_argsort = difference.argsort()
    difference_argsort = difference_argsort[~np.isnan(difference[difference_argsort])]
    difference_argsort = difference_argsort[::-1]
    difference_argsort = difference_argsort[:num_edges]
    difference = difference[difference_argsort]
    ct_edges = ct_edges[difference_argsort]

    # Format to df
    df = pd.DataFrame({'CT Edge': ct_edges, 'Difference': difference})

    # Plot
    sns.barplot(data=df, x='CT Edge', y='Difference', color='gray', ax=ax)
    plt.sca(ax)
    plt.xticks(rotation=90)
    plt.title(column)


def plot_ct_individual_edges(sid, *, ax, data, edges, heads, subject_ids, column=None, num_edges=15):
    # Defaults
    if column is None: column = get_attention_columns()[0]

    # Filter to subject
    data = data[:, :, np.array(subject_ids)==sid]

    # Filter to column
    data = data[:, np.array(heads)==column]

    # Filter to ct-ct edges
    ct_edges_mask = [string_is_synthetic(s.split(get_edge_string())[0]) and string_is_synthetic(s.split(get_edge_string())[1]) for s in edges]
    data = data[ct_edges_mask]
    ct_edges = np.array(edges)[ct_edges_mask]

    # Flatten and sort
    data = data.flatten()
    empty_mask = ~np.isnan(data)
    data = data[empty_mask]
    ct_edges = ct_edges[empty_mask]
    data_argsort = data.argsort()[::-1][:num_edges]

    # Filter
    data = data[data_argsort]
    ct_edges = ct_edges[data_argsort]

    # Format to df
    df = pd.DataFrame({'CT Edge': ct_edges, 'Attention': data})

    # Plot
    sns.barplot(data=df, x='CT Edge', y='Attention', color='gray', ax=ax)
    plt.sca(ax)
    plt.xticks(rotation=90)
    plt.title(column)


def plot_ct_individual_edge_comparison(sid, *, ax, data, edges, heads, subject_ids, columns=None, column_names=None, palette=['paleturquoise', 'salmon'], num_edges=15):
    # Defaults
    if columns is None: columns = get_attention_columns()[:2]

    # Filter to subject
    data = data[:, :, np.array(subject_ids)==sid].squeeze()

    # Filter to column
    data = data[:, [h in columns for h in heads]]

    # Filter to ct-ct edges
    ct_edges_mask = [string_is_synthetic(s.split(get_edge_string())[0]) and string_is_synthetic(s.split(get_edge_string())[1]) for s in edges]
    data = data[ct_edges_mask]
    ct_edges = np.array(edges)[ct_edges_mask]

    # Convert
    df = pd.DataFrame(data, index=ct_edges, columns=columns)
    df = df.reset_index(names='Edge').melt(id_vars='Edge', var_name='Head', value_name='Prioritization')
    if column_names is not None: df['Head'] = df['Head'].apply(lambda s: {c: cn for c, cn in zip(columns, column_names)}[s])

    # Filter
    mean_prioritization = df[['Edge', 'Prioritization']].groupby('Edge').sum().sort_values('Prioritization', na_position='first').loc[::-1]
    df = df.loc[df['Edge'].isin(mean_prioritization.index.to_numpy()[:num_edges])]

    # Plot
    sns.barplot(data=df, x='Edge', y='Prioritization', hue='Head', palette=['paleturquoise', 'salmon'], ax=ax)
    plt.sca(ax)
    plt.yscale('log')
    plt.xticks(rotation=90)


def plot_attention_histogram(sid, *, ax, data, edges, heads, subject_ids, column=None):
    # Defaults
    if column is None: column = get_attention_columns()[0]

    # Filter to subject
    data = data[:, :, np.array(subject_ids)==sid]

    # Filter to column
    data = data[:, np.array(heads)==column]

    # Flatten and filter
    data = data.flatten()
    empty_mask = ~np.isnan(data)
    data = data[empty_mask]
    new_edges = np.array(edges)[empty_mask]

    # Format to df
    df = pd.DataFrame({'Edge': new_edges, 'Attention': data})

    # Plot
    sns.histplot(data=df, x='Attention', color='gray', bins=30, kde=True, ax=ax)
    plt.sca(ax)
    plt.yscale('log')
    plt.xticks(rotation=90)
    plt.xlabel(column)


def create_subfigure_mosaic(shape_array, layout='constrained'):
    fig, axs = plt.subplots(shape_array.shape[0], shape_array.shape[1], layout=layout, figsize=(int((3/2) * shape_array.shape[1]), int((3/2) * shape_array.shape[0])))
    # Get gridspec
    gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    # Clear axes
    for a in axs.flatten(): a.remove()
    # Define subfigs
    subfigs = {}
    for key in np.unique(shape_array):
        # Skip blank spaces
        if key == '.': continue

        # Find bounding box (assume well-formed)
        tl = np.argwhere(shape_array==key).min(axis=0)
        br = np.argwhere(shape_array==key).max(axis=0) + 1

        # Record
        subfigs[key] = fig.add_subfigure(gridspec[tl[0]:br[0], tl[1]:br[1]])
    # Create subplots
    axs = {}
    for key, subfig in subfigs.items():
        axs[key] = subfig.add_subplot(1, 1, 1)
    # Formatting
    fig.set_constrained_layout_pads(w_pad=0, h_pad=0, wspace=.4, hspace=.4)  # *_pad is pad for figs (including subfigs), *_space is pad between subplots

    # NOTE: `bbox_inches='tight'` will eliminate whitespace from edges, even from '.' in mosaic.
    return fig, axs


def plot_ct_graph_from_sid(sid, *, ax, column=None, vertex_ids=None, g_pos=None, pos=None):
    # Defaults
    if column is None: column = get_attention_columns()[0]

    # Assemble
    g = compute_graph(load_graph_by_id(sid, column=column))

    # Filter to synthetic nodes
    g = filter_to_synthetic_vertices(g.copy(), vertex_ids=vertex_ids, depth=0)

    # Plot
    pos = get_graph_pos(g) if pos is None else convert_vertex_map(g_pos, g, pos)
    visualize_graph_base(
        g,
        pos=pos,
        ink_scale=.5,
        vertex_font_size=.1,
        edge_pen_width=gt.prop_to_size(g.ep.coef, .01, .05),
        mplfig=ax,
    )

    # Formatting
    ax.axis('off')

    # Return source graph and positions
    return {'g_pos': g, 'pos': pos}


def plot_prs_correlation(meta, *, data, edges, heads, subject_ids, ax=None, num_targets=10, min_samples=10, subsample=1., max_scale=False, df=None, prs_df=None, **kwargs):
    # Default
    if ax is None: ax = plt.gca()

    # Get data
    # TODO: Allow just passing prs_df
    if df is None: df = compute_prs_difference(meta, data=data, edges=edges, heads=heads, subject_ids=subject_ids, subsample=subsample, **kwargs)
    targets = (
        df[['edge', 'head']]
        .loc[df['n'] > min_samples]  # Filter by min_samples
        .iloc[:num_targets]  # Filter to top targets
        .to_numpy()
    )
    if prs_df is None: prs_df = get_prs_df(targets, meta=meta, data=data, edges=edges, heads=heads, subject_ids=subject_ids, **kwargs)

    # Max scale individual edges
    prs_df_format = prs_df.copy()
    if max_scale:
        for name in np.unique(prs_df_format['Name']):
            prs_df_format.loc[prs_df_format['Name']==name, 'Attention'] /= prs_df_format.loc[prs_df_format['Name']==name, 'Attention'].max()

    # Plot
    # NOTE: Right now, will not consider different edges on the same head
    prs_df_format = prs_df_format.rename(columns={'Risk': 'SCZ Risk'})
    sns.barplot(prs_df_format, x='Name', y='Attention', hue='SCZ Risk', hue_order=['low', 'mid', 'high'], ax=ax)
    # Formatting
    plt.sca(ax)
    plt.xticks(rotation=90)
    plt.xlabel(None)
    if max_scale:
        plt.ylabel('Prioritization (Max Scaled)')
    else:
        plt.ylabel('Prioritization')

    return df, prs_df


def plot_labels(axs, *, shape=None, **kwargs):
    """
    Plot labels for the given `axs`.  If `shape` is provided,
    alphabetize labels based on appearance order in `shape`.
    """
    if shape is not None: _, conversion, offset = alphabetize_shape(shape, return_dict=True, return_offset=True, **kwargs)
    for label, ax in axs.items():
        # NOTE: `in_layout=False` removes text from the constrained layout calculation
        ax.text(0, 1, conversion[label] if shape is not None else label, zorder=1, in_layout=False, ha='right', va='bottom', fontweight='bold', fontsize='xx-large', transform=ax.figure.transSubfigure)

    if shape is not None: return offset
