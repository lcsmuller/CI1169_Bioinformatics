import pandas as pd
from matplotlib import pyplot as plot
from sys import argv, exit

SEQUENCE_HITS_R5 = 'sequence_hits_r5.txt'
SEQUENCE_HITS_X4 = 'sequence_hits_x4.txt'

def generate_plots(path):
    # Load the datasets
    r5_data = pd.read_csv(f'{path}/{SEQUENCE_HITS_R5}', delim_whitespace=True, comment='#', header=None)
    x4_data = pd.read_csv(f'{path}/{SEQUENCE_HITS_X4}', delim_whitespace=True, comment='#', header=None)

    # Extract relevant columns for plotting
    r5_scores = r5_data[5]
    r5_evalues = r5_data[4]
    x4_scores = x4_data[5]
    x4_evalues = x4_data[4]

    # Convert evalues to numeric for plotting (they might be in scientific notation)
    r5_evalues = pd.to_numeric(r5_evalues, errors='coerce')
    x4_evalues = pd.to_numeric(x4_evalues, errors='coerce')

    # Scatter plot of E-values vs. Scores
    plot.figure(figsize=(12, 6))
    plot.scatter(r5_evalues, r5_scores, alpha=0.5, label='r5 Dataset', color='blue')
    plot.scatter(x4_evalues, x4_scores, alpha=0.5, label='x4 Dataset', color='red')
    plot.xscale('log')
    plot.xlabel('E-value')
    plot.ylabel('Score')
    plot.title('Scatter Plot of E-values vs. Scores')
    plot.legend()
    plot.text(0.95, 1.05, f'Path: {path}', transform=plot.gca().transAxes, ha='right', va='top')
    plot.savefig(f'{path}/scatter_plot.png')
    plot.close()

    # Box plot of Scores
    plot.figure(figsize=(10, 6))
    plot.boxplot([r5_scores, x4_scores], labels=['r5 Dataset', 'x4 Dataset'])
    plot.ylabel('Score')
    plot.title('Box Plot of Scores')
    plot.text(0.95, 1.05, f'Path: {path}', transform=plot.gca().transAxes, ha='right', va='top')
    plot.savefig(f'{path}/box_plot.png')
    plot.close()

if __name__ == "__main__":
    if len(argv) < 2:
        print("Please provide the path as a command line argument.")
        exit(1)

    generate_plots(argv[1])
