# How-to

1. Install dependencies
    ```sh
    sudo apt install hmmer clustalo
    ```
2. Align the sequences
    ```sh
    clustalo -i ccr5.fasta -o ccr5_aligned.fasta --force
    ```
3. Build HMM Profiles
    ```sh
    hmmbuild ccr5.hmm ccr5_aligned.fasta
    ```
4. Search for matches
    ```sh
    hmmsearch ccr5.hmm target_sequences.fasta > ccr5_search_results.txt
    ```
