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
    min_size = .7*g.num_vertices()**(-5/6)
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


def plot_legend(hub=False, ax=None):
    if not ax: ax = plt.gca()

    # Custom Legend
    # NOTE: Unfortunately, shapes need to be changed manually
    legend_elements = [
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='h', markerfacecolor=get_node_appearance('celltype')[0], label='Cell Type'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=get_node_appearance('tftg')[0], label='TF+TG'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='^', markerfacecolor=get_node_appearance('tf')[0], label='TF'),
        Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='o', markerfacecolor=get_node_appearance('tg')[0], label='TG'),
    ]
    if hub:
        legend_elements = [Line2D([0], [0], color='gray', linestyle='None', markersize=10, marker='8', markerfacecolor=get_node_appearance('hub')[0], label='Hub'),] + legend_elements
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


def plot_individual_edge_comparison(g, sample_ids, suffix='Attention Weights', highlight_outliers=True, ax=None):
    "Take concatenated graph `g` and plot a comparison between the original weights"
    if not ax: ax = plt.gca()

    # Assemble dataframe
    df = {k: [] for k in ['id'] + sample_ids}
    for e in g.edges():
        coefs = g.ep.coefs[e]
        df['id'].append(get_edge_string(g, e))
        for i, sample_id in enumerate(sample_ids):
            df[sample_id].append(coefs[i])
    df = pd.DataFrame(df)

    # Highlight outliers
    if highlight_outliers:
        # Get outliers
        difference = (df[sample_ids[1]] - df[sample_ids[0]]).to_numpy()
        outlier_idx, outlier_mask = get_outlier_idx(difference, return_mask=True)
        # Record
        df['outlier'] = outlier_mask
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
            ann = plt.text(
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
        # Make sure text doesn't overlap points or other text
        # NOTE: Doesn't work with plt.text right now
        # TODO: Tune for clearer visualization, perhaps figure out why arrows don't work

    # Plot
    sns.scatterplot(
        data=df,
        x=sample_ids[0],
        y=sample_ids[1],
        hue='outlier' if highlight_outliers else None,
        palette=['black', 'red'],
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

    # Plot y=x
    lims = [
        max(ax.get_xlim()[0], ax.get_ylim()[0]),
        min(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, '-', color='black', alpha=0.3)

    # Adjust text positions to avoid overlap
    # NOTE: Must be done after axis limit/scale changes
    if highlight_outliers:
        adjust_text(
            annotations,
            df[sample_ids[0]].to_numpy(),
            df[sample_ids[1]].to_numpy(),
            arrowprops=dict(arrowstyle='-', color='black', lw=.01),
        )


    # Broken
    # # Calculate discontinuities (From JAMIE)
    # max_dist = .2; pad = 1e-3
    # bounds = []
    # for vals in [df[sample_ids[0]], df[sample_ids[1]]]:
    #     bounds.append([])

    #     sorted_vals = np.sort(vals)
    #     min_val = sorted_vals[0]
    #     max_val = sorted_vals[0]
    #     for val in sorted_vals[1:]:
    #         if val - max_val > max_dist:
    #             bounds[-1].append((min_val - pad, max_val + pad))
    #             min_val = max_val = val
    #         else:
    #             max_val = val
    #     bounds[-1].append((min_val - pad, max_val + pad))

    # # Make broken plot
    # bax = brokenaxes(
    #     xlims=bounds[0],
    #     ylims=bounds[1],
    #     hspace=.15,
    #     wspace=.15,
    # )

    # # Get y=x
    # lims = [
    #     max(bounds[0][0][0], bounds[1][0][0]),
    #     min(bounds[0][-1][1], bounds[1][-1][1])]

    # # Plot everything
    # bax.plot(lims, lims, '-', color='black', alpha=0.3)
    # bax.scatter(df[sample_ids[0]], df[sample_ids[1]], color='black', edgecolors='white')

    return df


def get_mosaic(mosaic, scale=3):
    fig = plt.figure(figsize=(scale*len(mosaic[0]), scale*len(mosaic)), constrained_layout=True)
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
        ax.set_title(sid)
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
    for i, label in enumerate(lp.get_xticklabels()):
        if i % int(len(df)/num_x_labels) != 0: label.set_visible(False)


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

    # Plot
    # TODO: Set min size to represent true 0, get constant scaling
    scp = sns.scatterplot(
        data=df,
        x=column_name,
        y=index_name,
        hue=sign_name if multicolor else None,
        size=value_name,
        sizes=(0, size_max),
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
    # ax.axis('equal')
    ax.set_aspect('equal', 'box')
    # Zoom X
    margin = .5
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
