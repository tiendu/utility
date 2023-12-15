import argparse
import subprocess
import os
import tkinter as tk
from tkinter import filedialog, scrolledtext
import tkinter.simpledialog as simpledialog
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Entrez
import plotly.graph_objects as go
import threading
import webbrowser
import requests
from dataclasses import dataclass, field

# Ensure required packages are installed
def install_required_packages():
    required_packages = ["biopython", "plotly", "requests"]
    installed = subprocess.check_output(['pip', 'freeze']).decode().splitlines()
    installed_packages = [pkg.split('==')[0] for pkg in installed]
    for package in required_packages:
        if package not in installed_packages:
            subprocess.check_call(["pip", "install", package])
install_required_packages()

@dataclass
class AB1Viewer:
    root: tk.Tk
    file_path_var: tk.StringVar = field(default_factory=tk.StringVar)
    file_button: tk.Button = field(init=False)
    status_label: tk.Label = field(init=False)
    metadata_display: tk.Text = field(init=False)
    plot_button: tk.Button = field(init=False)
    blast_button: tk.Button = field(init=False)
    compare_button: tk.Button = field(init=False)
    alignment_display: scrolledtext.ScrolledText = field(init=False)
    results_display: scrolledtext.ScrolledText = field(init=False)
    blast_results: list = field(default_factory=list)

    def __post_init__(self):
        self.root.title("AB1 File Viewer")

        self.file_button = tk.Button(self.root, text="Select .ab1 File", command=self.load_file)
        self.file_button.pack(pady=10)

        self.status_label = tk.Label(self.root, text="", fg="blue")
        self.status_label.pack(pady=10)

        self.metadata_display = scrolledtext.ScrolledText(self.root, height=6, width=60)
        self.metadata_display.pack(pady=10)

        self.plot_button = tk.Button(self.root, text="Plot Chromatogram", command=self.plot_chromatogram)
        self.blast_button = tk.Button(self.root, text="Run BLAST Search", command=self.blast_and_display)
        self.plot_button.pack(pady=10)
        self.blast_button.pack(pady=10)

        self.results_display = scrolledtext.ScrolledText(self.root, height=20, width=80)
        self.results_display.pack(pady=10)

        self.compare_button = tk.Button(self.root, text="Compare with Best Hit", command=self.compare_with_best_hit)
        self.compare_button.pack(pady=10)

        self.alignment_display = scrolledtext.ScrolledText(self.root, height=20, width=80)
        self.alignment_display.pack(pady=10)

        self.help_button = tk.Button(self.root, text="Help", command=self.show_usage)
        self.help_button.pack(pady=10)

    def show_usage(self):
        usage = """
        AB1 File Viewer - Usage Guide

        1. Click 'Select .ab1 File' to choose an .ab1 file.
        2. Once selected, some metadata of the file will be displayed.
        3. Click 'Plot Chromatogram' to view the chromatogram of the selected file in a browser.
        4. Click 'Run BLAST Search' to get the BLAST results for the sequence in the selected file.
           The results, along with some annotations, will be displayed in the application.

        Note: Ensure you have an active internet connection when using the 'Run BLAST Search' feature.
        """
        # Create a new window for the usage
        usage_win = tk.Toplevel(self.root)
        usage_win.title("Usage Guide")
        usage_text = scrolledtext.ScrolledText(usage_win, height=20, width=80, wrap=tk.WORD)
        usage_text.insert(tk.END, usage)
        usage_text.pack(padx=10, pady=10)
        usage_text.config(state=tk.DISABLED)  # Make it read-only

    def load_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("AB1 files", "*.ab1")])
        if file_path:
            self.file_path_var.set(file_path)
            self.show_metadata(file_path)

    def show_metadata(self, file_name):
        record = SeqIO.read(file_name, "abi")
        annotations = record.annotations
        metadata = f"Sample Name: {annotations.get('sample_name', 'N/A')}\n"
        metadata += f"Machine Model: {annotations.get('machine_model', 'N/A')}\n"
        metadata += f"Run Start: {annotations.get('run_start', 'N/A')}\n"
        self.metadata_display.insert(tk.END, metadata)

    def plot_chromatogram(self):
        file_path = self.file_path_var.get()
        if file_path:
            html_path = create_plot(file_path)
            webbrowser.open('file://' + os.path.realpath(html_path))

    def blast_and_display(self):
        self.status_label.config(text="Running BLAST...")  # Indicate BLAST is running
        self.root.update_idletasks()  # Update the GUI

        file_path = self.file_path_var.get()
        record = SeqIO.read(file_path, "abi")
        sequence = str(record.seq)
        hits = blast_sequence(sequence)
        self.blast_results = hits   # Store the results

        for hit in hits:
            try:
                accession, title, length, e_value, similarity = hit
                annotations = get_uniprot_annotations(accession)

                # Creating bullet points for annotations
                bullet_content = ""
                if annotations:
                    bullet_content += f"• Entry Type: {annotations.get('entryType', 'N/A')}\n"
                    protein_desc = annotations.get('proteinDescription', {}).get('submissionNames', [{}])[0].get('fullName', {}).get('value', 'N/A')
                    bullet_content += f"• Protein Name: {protein_desc}\n"
                    bullet_content += f"• Primary Accession: {annotations.get('primaryAccession', 'N/A')}\n"
                    # You can continue adding fields here in similar fashion

                result_text = f"Accession: {accession}\nTitle: {title}\nLength: {length}\nE-value: {e_value}\nSimilarity: {similarity:.2f}%\n"
                result_text += "Annotations:\n" + bullet_content + "\n"

                self.results_display.insert(tk.END, result_text)
            except Exception as e:
                self.results_display.insert(tk.END, f"An error occurred while processing hit with accession {accession}: {str(e)}\n\n")
            self.status_label.config(text="")  # Clear the status label

    def compare_with_best_hit(self):
        self.status_label.config(text="Getting alignment...")  # Indicate alignment is running
        self.root.update_idletasks()  # Update the GUI

        email = simpledialog.askstring("Input", "Please enter your email for NCBI's API:")
        if not email:
            self.alignment_display.insert(tk.END, "No email provided. Cannot continue with BLAST.")
            return

        file_path = self.file_path_var.get()
        record = SeqIO.read(file_path, "abi")
        sequence = str(record.seq)

        if not self.blast_results:  # If no stored results, run BLAST
            hits = blast_sequence(sequence)
            self.blast_results = hits
        else:
            hits = self.blast_results  # Use stored results

        if not hits:  # No BLAST hits
            self.alignment_display.insert(tk.END, "No BLAST hits found to compare with.")
            return

        best_hit_accession, _, _, _, _ = hits[0]
        best_hit_sequence = fetch_sequence_by_accession(best_hit_accession, email)  # Provide the email to the function

        alignment = align_sequences(sequence, best_hit_sequence)
        self.alignment_display.insert(tk.END, alignment)
        self.status_label.config(text="")  # Clear the status label

def create_plot(file_name):
    record = SeqIO.read(file_name, "abi")
    traces = record.annotations['abif_raw']
    g_trace = traces['DATA9']
    a_trace = traces['DATA10']
    t_trace = traces['DATA11']
    c_trace = traces['DATA12']

    # Determine the number of nucleotides
    num_nucleotides = len(record.seq)
    # Set a pixel width per nucleotide
    pixel_width_per_nucleotide = 15
    # Calculate the total plot width
    plot_width = num_nucleotides * pixel_width_per_nucleotide

    fig = go.Figure()
    fig.add_trace(go.Scatter(y=g_trace, mode='lines', name='G', line=dict(color='green')))
    fig.add_trace(go.Scatter(y=a_trace, mode='lines', name='A', line=dict(color='red')))
    fig.add_trace(go.Scatter(y=t_trace, mode='lines', name='T', line=dict(color='blue')))
    fig.add_trace(go.Scatter(y=c_trace, mode='lines', name='C', line=dict(color='yellow')))

    # Determine an approximate threshold value
    max_intensity = max(max(g_trace), max(a_trace), max(t_trace), max(c_trace))
    threshold_value = 0.7 * max_intensity  # e.g., 70% of the max intensity

    fig.update_layout(
        title='AB1 Trace Data',
        xaxis_title='Position',
        yaxis_title='Signal Intensity',
        width=plot_width,
        legend=dict(x=-0.1, y=0.5),
        shapes=[
            dict(
                type='line',
                y0=threshold_value,
                y1=threshold_value,
                x0=0,
                x1=plot_width,
                line=dict(
                    color="Black",
                    width=2,
                    dash="dot"
                )
            )
        ]
    )

    if not os.path.exists("plots"):
        os.mkdir("plots")

    file_path = os.path.join("plots", "temp_plot.html")
    fig.write_html(file_path)
    return file_path

def blast_sequence(sequence, database="nr", max_hits=5):
    result_handle = NCBIWWW.qblast("blastn", database, sequence)
    blast_records = NCBIXML.parse(result_handle)

    hits = []
    for record in blast_records:
        for alignment in record.alignments[:max_hits]:
            # Consider the first HSP for simplicity (since BLAST results often include only one HSP per alignment)
            hsp = alignment.hsps[0]
            similarity_percentage = (hsp.identities / hsp.align_length) * 100
            hits.append((alignment.accession, alignment.title, alignment.length, hsp.expect, similarity_percentage))

    return hits

def get_uniprot_annotations(query, max_results=1):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    fields = "id,gene_names,protein_name,ec,go,cc_function,cc_pathway"

    params = {
        "query": query,
        "size": max_results,
        "fields": fields
    }

    headers = {
        "accept": "application/json"
    }

    response = requests.get(base_url, params=params, headers=headers)
    if response.status_code == 200:
        data = response.json()
        if "results" in data and len(data["results"]) > 0:
            return data["results"][0]  # Return first result
        else:
            return "No results found."
    else:
        return f"Error {response.status_code}: {response.text}"

def align_sequences(seq1, seq2):
    alignments = pairwise2.align.localxx(seq1, seq2)
    # For simplicity, return the first (best) alignment
    return format_alignment(*alignments[0])

def fetch_sequence_by_accession(accession, email):
    Entrez.email = email  # Use the provided email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)

if __name__ == "__main__":
    root = tk.Tk()
    viewer = AB1Viewer(root)
    root.mainloop()
