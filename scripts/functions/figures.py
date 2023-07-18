import matplotlib.pyplot as plt

from .plotting import *


### Figure functions
def figure_regulon_statistics(graph_summary):
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))

    # Phenotype and Sub-Phenotype vs TF Outgoing
    plot_tf_outgoing(graph_summary, ax=axs[0])
    # Phenotype and Sub-Phenotype vs Closeness
    plot_tf_closeness(graph_summary, ax=axs[1])

    return fig


def figure_diffusion(
        diff_g_individual_diffusion,
        diff_g_individual,
        diff_g_other,
        diff_graph_summary_coex,
        diff_graph_summary_att,
        *,
        meta,
        individual_sample_id,
        other_sample_id,
):
    # Generate figure layout
    scale = 3
    mosaic = [
        ['A1', 'A1', 'A1', 'A2', 'A2', 'A2', 'B1', 'B1', 'B1',],
        ['A1', 'A1', 'A1', 'A2', 'A2', 'A2', 'B1', 'B1', 'B1',],
        ['A1', 'A1', 'A1', 'A2', 'A2', 'A2', 'B1', 'B1', 'B1',],
        ['C1', 'C1', 'C1', 'C1', 'D1', 'D1', 'E1', 'E1', 'E1',],
        ['C1', 'C1', 'C1', 'C1', 'D1', 'D1', 'E2', 'E2', 'E2',],
    ]
    fig = plt.figure(figsize=(scale*len(mosaic[0]), scale*len(mosaic)), constrained_layout=True)
    axs = fig.subplot_mosaic(mosaic)
    axs['E1'].get_shared_x_axes().join(axs['E1'], axs['E2'])

    # Diffusion graph
    ax = axs['A1']
    plot_label('A', ax=ax)
    # ax.set_title('A', fontsize=30, weight='bold', loc='left')
    visualize_graph_diffusion(diff_g_individual_diffusion, ax=ax)
    ax.axis('off')
    ax.set_title('Diffusion')

    # Novel personal graph
    ax = axs['A2']
    visualize_graph(diff_g_individual, ax=ax)
    ax.axis('off')
    ax.set_title('Personal Subgraph')

    # Other sample personal graph
    ax = axs['B1']
    plot_label('B', ax=ax)
    visualize_graph(diff_g_other, ax=ax)
    ax.axis('off')
    ax.get_legend().set_visible(False)
    ax.set_title('Alternative Personal Subgraph')

    # Enrichment
    ax = axs['C1']
    plot_label('C', ax=ax)
    ax.axis('off')
    ax.set_title('Novel Gene Enrichment')
    # leave blank

    # NPS comparison
    ax = axs['D1']
    plot_label('D', ax=ax)
    plot_nps(meta, [individual_sample_id, other_sample_id], ax=ax)
    ax.set_title('NPS Comparison')

    # TF closeness coex
    ax = axs['E1']
    plot_label('E', ax=ax)
    plot_tf_closeness(diff_graph_summary_coex, ax=ax)
    ax.set_title('Average TF Closeness')
    ax.set_ylabel('Coexpression')
    ax.set_xlabel(None)
    ax.set_xticklabels([])

    # TF closeness att
    ax = axs['E2']
    plot_tf_closeness(diff_graph_summary_att, ax=ax)
    ax.get_legend().set_visible(False)
    ax.set_ylabel('Attention')

    return fig


def figure_data_driven(
        data_g_individual,
        diff_graph_summary,
        data_graph_summary,
        data_g_other,
        data_g_group,
):
    # Generate figure layout
    scale = 5
    mosaic = [
        ['A1', 'A1', 'B1', 'B1',],
        ['A1', 'A1', 'B2', 'B2',],
        ['C1', 'C1', 'D1', 'D1',],
        ['C1', 'C1', 'D1', 'D1',],
    ]
    fig = plt.figure(figsize=(scale*len(mosaic[0]), scale*len(mosaic)), constrained_layout=True)
    axs = fig.subplot_mosaic(mosaic)
    axs['B1'].get_shared_x_axes().join(axs['B1'], axs['B2'])

    # Personal graph
    ax = axs['A1']
    plot_label('A', ax=ax)
    visualize_graph(data_g_individual, ax=ax)
    ax.axis('off')
    ax.set_title('Personal Subgraph')

    # TF closeness diff
    ax = axs['B1']
    plot_label('B', ax=ax)
    plot_tf_closeness(diff_graph_summary, ax=ax)
    ax.set_ylabel('Diffusion')
    ax.set_xlabel(None)
    ax.set_xticklabels([])

    # TF closeness data
    ax = axs['B2']
    plot_tf_closeness(data_graph_summary, ax=ax)
    ax.set_title(None)
    ax.get_legend().set_visible(False)
    ax.set_ylabel('Data Driven')

    # Personal graph
    ax = axs['C1']
    plot_label('C', ax=ax)
    visualize_graph(data_g_other, ax=ax)
    ax.axis('off')
    ax.get_legend().set_visible(False)
    ax.set_title('Alternate Personal Subgraph')

    # Visualization of novel subnetworks beyond diffused
    ax = axs['D1']
    plot_label('D', ax=ax)
    visualize_graph(data_g_group, ax=ax)
    # visualize_graph_state(data_g_group, ax=ax)
    ax.axis('off')
    ax.set_title('Group Subgraph')

    return fig
