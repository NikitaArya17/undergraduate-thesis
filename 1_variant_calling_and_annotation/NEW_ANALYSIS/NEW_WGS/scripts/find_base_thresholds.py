import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Analyze VCF metrics for GATK Bootstrapping.')
    parser.add_argument('-i', '--input', required=True, help='Input table from VariantsToTable')
    parser.add_argument('-t', '--type', required=True, choices=['snps', 'indels'], help='Mutation type')
    parser.add_argument('-p', '--percentile', type=float, default=0.90, help='Percentile for "Gold Set" (default 0.90)')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    
    metrics = ['QD', 'MQ', 'FS', 'SOR', 'MQRankSum', 'ReadPosRankSum']
    
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f'Distribution of Quality Metrics: {args.type.upper()}', fontsize=20)
    axes = axes.flatten()

    results = []

    for i, col in enumerate(metrics):
        if col not in df.columns:
            continue
            
        data = df[col].dropna()
        
        if col in ['QD', 'MQ']:
            thresh = data.quantile(args.percentile)
        else:
            # For all other metrics, lower values are higher quality
            thresh = data.quantile(1 - args.percentile)
        
        results.append(f"{col}: {thresh:.4f}")

        sns.histplot(data, kde=True, ax=axes[i], color='teal' if args.type == 'snps' else 'indianred')
        axes[i].axvline(thresh, color='red', linestyle='--', label=f'Threshold ({thresh:.2f})')
        axes[i].set_title(f'Distribution of {col}')
        axes[i].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"{args.type}_metrics_dist.png")
    
    print(f"--- {args.type.upper()} BOOTSTRAP THRESHOLDS ({args.percentile*100}th Percentile) ---")
    print("\n".join(results))
    print(f"--- Plot saved as {args.type}_metrics_dist.png ---")

if __name__ == "__main__":
    main()