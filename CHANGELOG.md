### 2023-10-17
- Adjust node sizes
- Adjust heatmap size
- Temporarily adjust heatmap scale

### 2023-10-16
- Clarified comparisons
- Miscellaneous bugfixes
- Modified aspect ratios
- New data
- New individual dosage analysis
- New visualizations for figures 3 and 4
- Other small updates

### 2023-10-13
- Figure revisions

### 2023-10-11
- Figure revisions

### 2023-09-19 (1-2)
- Additional plot input sanitization
- Figure revisions
- Fix enrichment plot scaling to start from zero
- New FREEZE and graphs with corresponding reading adjustments
- Recalculation of manual notebook parameters

### 2023-09-18 (1-2)
- Add FDR correction to SNP analysis
- Add eQTL SNP visualization to individual plots
- Calculate and read real enrichments
- Figure revisions
- New pathway enrichment for individual comparisons
- Refine SNP analysis and add utility arguments
- Swap enrichment axes

### 2023-09-15
- Add multiple Manhattan plots

### 2023-09-13
- Add various utility functions
- Implement dosage analysis

### 2023-09-12 (1-2)
- Add gene prioritization for group analyses
- Add framework for dosage correlation analysis
- Add stacked barchart visualization for `plot_prediction_confusion`
- Figure revisions
- Sanitize input and optimize processing for `compute_prediction_confusion`
- Separate combined figures into individual elements
- Small logic and formatting changes

### 2023-09-06 (1-2)
- Add legend titles
- Additional formatting for `plot_BRAAK_comparison`
- Additional utility plotting functions for legend, including `plot_remove_legend` and `plot_outside_legend`
- Consistent coloration in figure 3
- Figure updates
- Fix bug in new node attribute calculation function
- Fix zoom in `plot_individual_edge_comparison`
- Fixed text coloration in graph plotting
- Highlight outliers in `plot_individual_edge_comparison`
- More organization for notebook
- Several small formatting fixes

### 2023-09-05 (1-4)
- Add `get_node_attributes` function to sync color/shape/size changes between plots
- Figure layout and plot format updates
- Legend recalculation
- Update graph node colors to match color dictionary

### 2023-08-31 (1-2)
- Adapt `plot_contrast_curve` to additional DataFrame-type arguments
- Add batch graph loading with `get_graphs_from_sids`
- Additional arguments for many functions
- Additional population analysis plots
- Individual figure updates
- SNP analysis plot framework and dummy figure

### 2023-08-30 (1-2)
- Add general-purpose circle-style heatmap plotting with `plot_circle_heatmap`
- Add head comparison plot, `plot_head_comparison`
- Add many utility functions
- Add multi-column capability to `load_graph_by_id`
- Add supplemental figures
- Add utility function `get_attention_columns`
- Change default TF-TG recalculation timing in `concatenate_graphs` to `after`, also add configurable option
- Change filename structure
- Change plot ordering and grouping
- Change several plot styles
- Figure updates
- Fix several small bugs and unoptimal behaviors
- Flip enrichment plot
- Implement `mean` sorting method for `plot_contrast_curve`
- Manual plot adjustments
- Optimize `remove_duplicate_edges` by using hash table for duplicate edges (I am dumb)
- Revise automatic threshold for `concatenate_graphs`
- Revised individual loading structure for `graph_analysis` notebook
- Small plot revisions for readability

### 2023-08-29 (1-3)
- Add more objectives for figure 4 plots
- Added filtering for aggregate graphs
- Fixed `concatenate_graphs` bug with `threshold=True` filtering out all edges
- Increase flexibility for many (previously) BRAAK calculations in figure 4
- New figures and layouts
- New plot computations
- Place new plots into functions
- Various bugfixes and small formatting optimizations

### 2023-08-28
- Add many utility functions
- Implement new figure plots
- Label plots in notebook
- Many more arguments for a variety of functions
- Many alternative new plots and layouts
- `concatenate_graphs` bugfixes and optimizations

### 2023-08-23
- Figure updates

### 2023-08-21 (1-3)
- Add features and automatic culling to `compute_aggregate_edge_summary` and `compute_edge_summary`
- Add function to filter graph based on synthetic vertices
- Add individual plot visualizations for specific subjects
- Add new plots and computations for characterization of contrast subgroups
- Additional plotting functionality for `plot_individual_edge_comparison`, allowing for broken axes
- Fix individual edge comparison plot bug causing common edges to be excluded
- Update edge text translation

### 2023-08-17
- Additional automated filtering capabilities
- New scaling strategy, should be more consistent

### 2023-08-16
- Add common edge detection, and set as default for concatenation
- Change isolated node visualization for `plot_graph_comparison`
- Everything generally looks better
- Fix scaling issues with position calculations
- Fix vertex attribute updating upon concatenation
- General bugfixes
- Make graph computation more modular
- Make node sizes vary by type
- New custom plot scaling methods to counter auto-scaling by `graph-tool`
- New synthetic node detection

### 2023-08-15 (3)
- Add new protocol for `cull_isolated_leaves`
- Fix for gene prioritization cell gene isolation
- Fix `num_x_labels` bug
- Fix `remove_text_by_centrality` `efilt` bug
- Full runs for aggregate comparisons
- Greatly optimize `concatenate_graphs` through `_remove_duplicate_edges`
- Many small optimizations and visualizations
- Optimize `compute_edge_summary`
- Optimize `plot_aggregate_edge_summary`
- Visibility updates

### 2023-08-15 (1-2)
- Added function to compare individual and aggregate graphs
- Move figure code to notebooks for the future
- New graph comparison plots
- Visibility updates

### 2023-08-08 (1-2)
- Changed the majority of plots
- Many new plots, including enrichment and individual edge comparison
- New strategy for plot exporting
- Refactored computation section
- Update figures

### 2023-08-01 (1-3)
- Add coloration option to graph combination
- Add example plots
- Figure alternatives
- Increase compatibility of `concatenate_graphs` with visualization and utility functions
- Sankey update

### 2023-07-25
- Add graph embedding analysis

### 2023-07-24 (1-2)
- Add enhanced graph subsetting features
- Add filtering to graph computation
- Add graph embedding loading
- Added several functions facilitating common gene positioning across graph plots
- Adjust coloration and formatting
- Fix PDF line weights by using alpha
- Fix certain plot labels, especially for NPS comparison
- Make file loading more modular
- Revise figures
- Scale text properly for `visualize_graph_diffusion`
- Tuning of variable selection for data figure
- Utilize new data

### 2023-07-18
- Change file structure

### 2023-07-13 (1-5)
- Change to PDF
- First figure versions for `diffusion` and `data`

### 2023-07-12
- Reorganization
- New plot format

### 2023-07-10
- Add new plots
- Code refactor in `General Analysis`
- Many new functions for plotting and computation (e.g. `subset_graph`)
- Reformat old plots
- Switch to mosaic plotting

### < 2023-07-10
- Creation of the changelog
