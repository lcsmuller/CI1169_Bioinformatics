import pandas as pd

# Define the path to the GTF file
gtf_file = 'gencode.v46.annotation.gtf'

# Columns of interest
columns = ['gene_id', 'gene_symbol', 'chromosome', 'start', 'end']

# Create an empty DataFrame
gene_positions = pd.DataFrame(columns=columns)

# Read the GTF file and extract gene positions
with open(gtf_file, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] == 'gene':
            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            info = {entry.split()[0]: entry.split()[1].strip('"') for entry in fields[8].split(';') if entry}
            gene_id = info.get('gene_id')
            gene_name = info.get('gene_name')
            if gene_id and gene_name:
                gene_positions = gene_positions._append({
                    'gene_id': gene_id,
                    'gene_symbol': gene_name,
                    'chromosome': chrom,
                    'start': start,
                    'end': end
                }, ignore_index=True)

# Save to CSV
gene_positions.to_csv('gene_positions.csv', index=False)
