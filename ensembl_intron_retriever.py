import argparse
import requests
from concurrent.futures import ProcessPoolExecutor

def fetch_transcript_info(gene_name):
    # Define the Ensembl REST API server URL
    server = "https://rest.ensembl.org"

    # Build the API endpoint to look up gene information using the species and gene name
    gene_endpoint = f"{server}/lookup/{gene_name}?expand=1"

    # Make the API request to get gene information
    response = requests.get(gene_endpoint, headers={"Content-Type": "application/json"})

    # Check if the request was successful
    if not response.ok:
        print(f"Error: {response.text}")
        return None

    # Parse the JSON response containing gene data
    gene_data = response.json()

    # Get the gene name
    gene_name = gene_data.get("display_name")

    # Get the species information
    species = gene_data.get("species")

    # Get the list of transcripts for the gene
    transcripts = gene_data.get("Transcript")

    # Check if any transcripts are found for the gene
    if not transcripts:
        print(f"No transcripts found for the gene {gene_name} ({species}).")
        return None

    # Process intron positions for each transcript
    results = []
    for transcript in transcripts:
        # Get the transcript name/variant
        transcript_name = transcript.get("id")

        # Get the strand information
        strand = "+" if transcript.get("strand") == 1 else "-"

        # Get the list of exons for the transcript
        exons = transcript.get("Exon")

        # Check if any exons are found for the transcript
        if not exons:
            print(f"No exons found for the gene transcript {transcript_name} of {gene_name} ({species}).")
            continue

        # Calculate intron positions based on exon positions
        introns = []
        for i in range(len(exons) - 1):
            exon_end = exons[i]["end"]
            next_exon_start = exons[i + 1]["start"]
            intron_start = exon_end + 1
            intron_end = next_exon_start - 1
            introns.append({"start": intron_start, "end": intron_end})

        # Store the intron positions for the transcript
        results.append((gene_name, species, transcript_name, strand, introns))

    return results

def get_introns(gene_names, threads):
    # Use ProcessPoolExecutor to fetch transcript info concurrently
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_results = {executor.submit(fetch_transcript_info, gene_name): gene_name for gene_name in gene_names}

        # Process the completed futures and print intron positions
        for future in future_to_results:
            gene_name = future_to_results[future]
            results = future.result()

            if results:
                for gene_name, species, transcript_id, strand, introns in results:
                    header = f"Intron Positions for Gene: {gene_name} ({species}), Transcript: {transcript_id}, Strand: {strand}"
                    print(header)
                    for intron in introns:
                        print(f"Chromosome: {transcript_id}, Start: {intron['start']}, End: {intron['end']}")
                    print("=" * len(header))

def main():
    # Set up argparse to parse the command-line arguments
    parser = argparse.ArgumentParser(description="Retrieve intron positions for genes from Ensembl.")
    parser.add_argument('-i', '--input_genes', nargs='+', type=str, required=True,
                        help="Gene names or symbols separated by space (e.g., ENSG00000141510 ENSG00000145335)")
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help="Number of threads/cores to use for processing (default: 4)")
    args = parser.parse_args()

    # Call the function to get intron positions for the specified genes and species
    get_introns(args.input_genes, args.threads)
    
if __name__ == "__main__":
    main()
