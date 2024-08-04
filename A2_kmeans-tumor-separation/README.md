# How-to

1. Install dependencies
    ```sh
    # For manipulating single-cell data and clustering
    pip3 install infercnvpy scanpy scikit-learn scikit-misc
    ```
2. Download GTF file from GENCODE
    ```sh
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
    gunzip gencode.v46.annotation.gtf.gz
    ```
3. Run build
    ```sh
    python3 extract_gene_positions.py
    python3 generate_plots.py
    ```
