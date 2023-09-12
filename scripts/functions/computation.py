import itertools as it

import graph_tool.all as gt
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import SGDClassifier
from sklearn.neural_network import MLPClassifier
from sklearn import metrics
from sklearn.model_selection import train_test_split
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
    print('Calculating statistics...')
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
    # Filter to high-valued edges
    if filter:
        graph = graph[graph['coef'] >= graph['coef'].quantile(filter)]

    # Convert graph to graph_tools
    list_of_tuples = list(graph.itertuples(index=False, name=None))
    g = gt.Graph(list_of_tuples, hashed=True, eprops=[('coef', 'double')])

    # Label self loops and add color
    g.ep.self_loop = g.new_edge_property('bool')
    gt.label_self_loops(g, eprop=g.ep.self_loop, mark_only=True)
    g.vp.self_loop_value = g.new_vertex_property('double')
    g.ep.color = g.new_edge_property('vector<double>')
    for e in g.edges():
        # Label self-loops
        if g.ep.self_loop[e]:
            g.vp.self_loop_value[e.source()] = g.ep.coef[e]

        # Add color to edges
        alpha = get_alpha(g.ep.coef[e])
        g.ep.color[e] = [0, 0, 0, alpha]

    # Determine color and flavor text
    g = assign_vertex_properties(g)

    # View without self-loops
    g_noself = gt.GraphView(g, efilt=lambda e: not g.ep.self_loop[e])

    return g_noself


def compute_edge_summary(graphs=None, concatenated_graph=None, *, subject_ids, min_common_edges=1, threshold=None):
    "Return edge x subject dataframe for all graphs"
    # TODO: Make `subject_ids` not required
    # Setup
    assert graphs is not None or concatenated_graph is not None
    if graphs is not None: concatenated_graph = concatenate_graphs(*graphs, threshold=threshold)

    # Format edge weights
    # df = pd.DataFrame(columns=['Edge']+subject_ids)
    df = {k: [] for k in ['Edge']+subject_ids}
    print('Collecting edges...')
    for e in tqdm(concatenated_graph.edges(), total=concatenated_graph.num_edges()):
        edge_name = get_edge_string(concatenated_graph, e)
        coefs = concatenated_graph.ep.coefs[e]
        # Take only edges which are common between two or more graphs
        if sum([c!=0 for c in coefs]) >= min_common_edges:
            row = [edge_name] + list(coefs)
            # df.loc[df.shape[0]] = row  # Slow
            for k, v in zip(df, row):
                df[k].append(v)
    df = pd.DataFrame(df)

    # Find variance and mean
    df['Variance'] = np.var(df.iloc[:, 1:1+len(coefs)], axis=1)
    df['Mean'] = np.mean(df.iloc[:, 1:1+len(coefs)], axis=1)
    # df = df.sort_values(['Variance', 'Mean'], ascending=[True, False])
    # df['index'] = list(range(len(df)))  # Store sort

    return df, concatenated_graph


def compute_aggregate_edge_summary(contrast_subject_ids, *, column, max_graphs=np.inf, threshold=None):
    "Return concatenated graphs for a contrast"
    # TODO: Add list input for `contrast_subject_ids`
    # For each subgroup of the contrast
    contrast_concatenated_graphs = {}; contrast_concatenated_subject_ids = {}
    for key, subject_ids in contrast_subject_ids.items():  # , total=len(contrast_subject_ids)
        # Get concatenated graph
        graphs = []; sids = []
        for sid in np.random.choice(subject_ids, len(subject_ids), replace=False):  # , total=len(subject_ids)
            try: graphs.append(compute_graph(load_graph_by_id(sid, column=column)))
            except: continue
            sids.append(sid)
            if len(sids) >= max_graphs: break
        # Skip if no graphs found
        if len(sids) == 0:
            continue
        concatenated_graph = concatenate_graphs(*graphs, threshold=threshold)
        contrast_concatenated_graphs[key] = concatenated_graph
        contrast_concatenated_subject_ids[key] = sids

        # Cleanup
        del graphs

    return contrast_concatenated_graphs, contrast_concatenated_subject_ids


def compute_contrast_summary(contrast, *, column, population=True):
    "Return dataframe with edge mean and variance for all subgroups in contrast"
    # TODO: Allow for multiple columns at the same time, requires overhaul of multiple functions, including concat...
    # Get subgroup variance
    contrast_concatenated_graphs, contrast_concatenated_subject_ids = compute_aggregate_edge_summary(
        contrast_subject_ids=contrast, column=column)  # Threshold included for reliable variance calculation
    df_subgroup = {}
    for subgroup in contrast_concatenated_graphs:
        df, concatenated_graph = compute_edge_summary(
            concatenated_graph=contrast_concatenated_graphs[subgroup], subject_ids=contrast_concatenated_subject_ids[subgroup])
        df = df[['Edge', 'Mean', 'Variance']]
        df['Subgroup'] = subgroup
        df_subgroup[subgroup] = df

    # Get population variance
    if population:
        sample_ids = {'Population': sum([contrast[k] for k in contrast], [])}
        contrast_concatenated_graphs, contrast_concatenated_subject_ids = compute_aggregate_edge_summary(
            contrast_subject_ids=sample_ids, column=column)  # Threshold included here for filtering
        df, concatenated_graph = compute_edge_summary(
            concatenated_graph=contrast_concatenated_graphs['Population'], subject_ids=contrast_concatenated_subject_ids['Population'])
        df = df[['Edge', 'Mean', 'Variance']]
        df['Subgroup'] = 'Population'
        df_subgroup['Population'] = df  # Add population as subgroup

    return df_subgroup


def compute_BRAAK_comparison(
        contrast,
        *,
        meta,
        column,
        target='BRAAK_AD',
        edges_include=None,
        edge_percentile=90,
        num_edges=5,
        seed=42):
    """
    Compute df with attention scores of `num_edges` random edges with edge commonality in the
    `edge_percentile` percentile annotated by `target` in `meta` over all individuals in the
    contrast.
    """
    # Calculate
    sids = sum([sids for _, sids in contrast.items()], [])
    all_graphs, sids = load_many_graphs(sids, column=column)
    all_graphs = [compute_graph(graph) for graph in all_graphs]
    df, _ = compute_edge_summary(graphs=all_graphs, subject_ids=sids)

    # Process
    df = df.drop(columns=['Variance', 'Mean'])
    df = pd.melt(df, id_vars=['Edge'], var_name='Subject ID', value_name='Attention')
    df.index = df['Subject ID']
    df_meta = meta.copy()[[target]]
    df_meta.index = meta['SubID']
    df = df.join(df_meta, how='left').reset_index(drop=True)

    # Format
    df = df.loc[df['Attention'] != 0]  # Remove 0 attention
    if edges_include is None:
        all_possible_edges, counts = np.unique(df['Edge'], return_counts=True)
        all_possible_edges = all_possible_edges[counts > np.percentile(counts, edge_percentile)]
        # TODO: Use highest variance or similar rather than random edges
        np.random.seed(seed)
        edges_include = np.random.choice(all_possible_edges, num_edges, replace=False)
    df = df.loc[[e in edges_include for e in df['Edge']]]

    return df, edges_include


def compute_prediction_confusion(
        contrast,
        *,
        meta,
        column,
        prioritized_edges,
        target='BRAAK_AD',
        classifier_type='SGD',
        random_state=42):
    # Calculate
    sids = sum([sids for _, sids in contrast.items()], [])  # All sids in contrast
    all_graphs, sids = load_many_graphs(sids, column=column)
    all_graphs = [compute_graph(graph) for graph in all_graphs]
    df, concatenated_graph = compute_edge_summary(graphs=all_graphs, subject_ids=sids)

    # Filter
    df = df.drop(columns=['Variance', 'Mean'])
    df = df.loc[[e in prioritized_edges for e in df['Edge']]]

    # Format
    X = np.array(df)[:, 1:].T
    df_meta = meta.copy()
    df_meta.index = df_meta['SubID']
    df_meta = df_meta.loc[list(df.columns)[1:]].reset_index(drop=True)
    y = df_meta[target].to_numpy().astype(str)

    # Remove nan
    is_nan = pd.isna(df_meta[target])
    X, y = X[~is_nan], y[~is_nan]

    # Remove classes with too few samples
    # NOTE: This prevents an error with multiclass prediction for sklearn
    min_samples = 2
    y_unique, y_counts = np.unique(y, return_counts=True)
    indices_to_include = sum([y==val for val, count in zip(y_unique, y_counts) if count >= min_samples])
    indices_to_include = np.array(indices_to_include).astype(bool)
    X, y = X[indices_to_include], y[indices_to_include]

    # Formatting and object standardization
    names = np.unique(y)

    # Predict
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=random_state)
    if classifier_type == 'SGD':
        classifier = SGDClassifier(random_state=random_state).fit(X_train, y_train)
    elif classifier_type == 'MLP':
        classifier = MLPClassifier(random_state=random_state).fit(X_train, y_train)
    else:
        raise ValueError(f'Classifier type {classifier_type} not found.')
    accuracy = classifier.score(X_test, y_test)

    # Evaluate
    confusion_matrix = metrics.confusion_matrix(y_test, classifier.predict(X_test), labels=names)
    row_sum = confusion_matrix.sum(axis=1)
    row_acc = np.diag(confusion_matrix) / row_sum
    col_sum = confusion_matrix.sum(axis=0)
    col_acc = np.diag(confusion_matrix) / col_sum
    df = pd.DataFrame(
        confusion_matrix,
        # Predicted names
        columns=[
            f'Predicted {name} (n={int(n)}, acc={racc:.3f})'
            for name, n, racc in zip(names, col_sum, col_acc)],
        # True names
        index=[
            f'{name} (n={int(n)}, acc={racc:.3f})'
            for name, n, racc in zip(names, row_sum, row_acc)])

    return df, accuracy


def compute_head_comparison(subject_ids, **kwargs):
    "Return edge x head difference DataFrame"
    # Setup
    all_columns = get_attention_columns()

    # Get graphs
    joined_graphs = get_many_graph_lists(subject_ids, all_columns)

    # CLI
    print(f'{joined_graphs.shape[0]} common edges found')

    # Calculate differences
    for column in all_columns:
        joined_graphs[column] = joined_graphs[column+'_s1'] - joined_graphs[column+'_s2']

    # Get top idx
    idx_to_include = get_top_idx(joined_graphs.abs(), all_columns, **kwargs)

    # Filter
    joined_graphs = joined_graphs.iloc[idx_to_include][all_columns]

    return joined_graphs


def get_graphs_from_sids(subject_ids, *, method='attention', column=None):
    # TODO: Add try, except
    if method == 'coex':
        return [
            cull_isolated_leaves(
                compute_graph(
                    scale_edge_coefs_list(
                        load_graph_by_id(
                            sid,
                            source='coexpression'),
                        1./60),
                    filter=.9))
            for sid in subject_ids]
    elif method == 'attention':
        return [compute_graph(load_graph_by_id(sid, column=column)) for sid in subject_ids]
