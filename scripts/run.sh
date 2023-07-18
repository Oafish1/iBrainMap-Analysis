eval "$(conda shell.bash hook)"
conda activate GNN
jupyter nbconvert --execute --to notebook --inplace graph_analysis.ipynb
