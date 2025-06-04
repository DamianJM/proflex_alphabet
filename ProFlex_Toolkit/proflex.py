 # Global Imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os, csv, pickle
from prody import parsePDB, ANM, calcSqFlucts
from collections import defaultdict
from itertools import islice

from Bio.PDB import PDBParser, Polypeptide, Superimposer, CaPPBuilder
from Bio import Align, pairwise2
from Bio.pairwise2 import format_alignment



try:
    import pymol2
except ImportError:
    pymol2 = None
    print("Pymol2 is not available which may be a consequence of failed pip installation")
    print("You may need to install pymol via a specific wheel. Please consult the github repo for instructions.")

class ProFlex:
    # Proflex empirically defined percentiles from global binning
    PERCENTILES: list = [
        0.0, 0.0002358430342649322, 0.00043512998728183964, 0.0006312815003185825, 0.0008351154949256513,
        0.001051505902672258, 0.001284168879166294, 0.0015360183731002995, 0.0018095910346847563, 0.0021085739460822276,
        0.002435563233031076, 0.002794854968964461, 0.0031904817416206366, 0.003626731382937652, 0.00410782947837465,
        0.004639059751976868, 0.005227352083684521, 0.005875978579139607, 0.006592542233410673, 0.007385499844320906,
        0.008262005819282724, 0.009231685364363246, 0.010302705256334073, 0.01148654709889712, 0.012795053494303177,
        0.014238649728642557, 0.015830891978568143, 0.017587869426746856, 0.01952463798475797, 0.02166113057134577,
        0.024026338047061186, 0.02663622503402896, 0.029512027757647025, 0.03269764330871081, 0.03623155292367679,
        0.0401535951145026, 0.04452904962775647, 0.04944217442639428, 0.05497012728148079, 0.06121154832082872,
        0.06834379731723829, 0.07659002272378347, 0.08621519073471832, 0.09756543518478902, 0.11118477495059192,
        0.12791181857140152, 0.14895298232149787, 0.17623308052226272, 0.2137183434951446, 0.26921283921832057,
        0.36244452191822385, 0.5, 1.0
    ]

    ALPHABET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    SUBSTITUTION_MATRIX, ROWS, COLUMNS = [], [], []

    with open('proflex_substitution_matrix.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        COLUMNS = next(reader)[1:]
        for row in reader:
            ROWS.append(row[0])
            numeric_row = [float(cell) for cell in row[1:] if cell.strip() != '']
            SUBSTITUTION_MATRIX.append(numeric_row)

    SUBSTITUTION_MATRIX = np.array(SUBSTITUTION_MATRIX)

    def __init__(self):
        """Initialize the ProFlex class."""
        pass

    @staticmethod
    def scale_rmsf(values):
        """
        Scale RMSF values between 0 and 1.
        """
        min_value = min(values)
        max_value = max(values)
        scaled_values = [(val - min_value) / (max_value - min_value) for val in values]
        return scaled_values

    def translate_to_alphabet(self, scaled_values):
        """
        Translate scaled RMSF values to alphabet based on predefined percentiles.
        """
        n_bins = len(self.ALPHABET)
        bin_indices = np.digitize(scaled_values, self.PERCENTILES) - 1
        bin_indices = np.clip(bin_indices, 0, n_bins - 1)
        mapped_letters = [self.ALPHABET[idx] for idx in bin_indices]
        return mapped_letters

    def backtranslate_to_values(self, letters):
        """
        Backtranslate alphabet characters to scaled RMSF values based on percentiles.
        """
        n_bins = len(self.ALPHABET)
        bin_widths = np.diff(self.PERCENTILES)
        bin_midpoints = self.PERCENTILES[:-1] + bin_widths / 2
        letter_to_value = {self.ALPHABET[i]: bin_midpoints[i] for i in range(n_bins)}
        values = [letter_to_value[letter] for letter in letters]
        return values

    def encode_sequence(self, rmsf_values):
        """
        Encode RMSF vector to ProFlex
        """
        scaled_values = self.scale_rmsf(rmsf_values)
        return self.translate_to_alphabet(scaled_values)

    def decode_sequence(self, alphabet_sequence):   # note that percentile centre points are used as values in decoding
        """
        Derive scaled RMSF from ProFlex
        """
        return self.backtranslate_to_values(alphabet_sequence)

    @staticmethod
    def read_rmsf_fasta(file_path):
        """
        Reads FASTA and imports RMSF values
        """
        sequences = {}
        current_header = None

        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line[1:]  # Remove the '>' from the header
                    sequences[current_header] = []
                elif current_header is not None and line:
                    try:
                        # Convert comma-separated values to a list of floats
                        values = list(map(float, line.split(',')))
                        sequences[current_header].extend(values)
                    except ValueError as e:
                        print(f"Error parsing line under header {current_header}: {e}")

        return sequences

    def process_multifasta(self, seqs):
        """
        Process multiple FASTA-like sequences by encoding and decoding them.
        """
        for header, rmsf_values in seqs.items():
            encoded_sequence = self.encode_sequence(rmsf_values)
            print(f"{header} Encoded Sequence:", encoded_sequence)

            decoded_values = self.decode_sequence(encoded_sequence)
            print(f"{header} Decoded Values:", decoded_values)

    @staticmethod
    def nma_analysis_to_dict(pdb_file):
        """
        Perform Normal Mode Analysis (NMA) on a PDB file and return the RMSF values in a dictionary format.

        Parameters:
        pdb_file (str): Path to the PDB file.

        Returns:
        dict: A dictionary where the key is the PDB title and the value is a list of RMSF values.
        """
        # Parse the PDB file
        protein = parsePDB(pdb_file)

        # Select the alpha carbons (CA) for the ANM model (NMA usually focuses on backbone atoms)
        calphas = protein.select('calpha')

        # Perform ANM (Anisotropic Network Model) analysis
        anm = ANM('ANM analysis')
        anm.buildHessian(calphas)
        anm.calcModes(n_modes=20)  # Calculating the first 20 modes

        # Calculate the squared fluctuations from the modes
        sq_fluctuations = calcSqFlucts(anm[:20])  # Calculate fluctuations for the first 20 modes

        # Convert squared fluctuations to RMSF by taking the square root
        rmsf = np.sqrt(sq_fluctuations)

        # Create a dictionary with the PDB title as the key and the RMSF values as the value
        rmsf_dict = {protein.getTitle(): rmsf.tolist()}  # Convert RMSF values to a list

        return rmsf_dict
   
    def print_substitution_matrix(self):
        print(self.SUBSTITUTION_MATRIX)

    def visualise_substitution_matrix(self):
        plt.figure(figsize=(10, 8))
        plt.imshow(self.SUBSTITUTION_MATRIX, cmap='viridis', aspect='auto')
        plt.colorbar(label='Log Odds')
        plt.xticks(ticks=np.arange(len(self.COLUMNS)), labels=self.COLUMNS, rotation=90)
        plt.yticks(ticks=np.arange(len(self.ROWS)), labels=self.ROWS)
        plt.xlabel('')
        plt.ylabel('')
        plt.title('ProFlex Substitution Matrix')
        plt.tight_layout()
        plt.show()
   
class NGramDatabase:
    def __init__(self, n=5):
        self.ngram_database = defaultdict(list)
        self.sequence_map = {}
        self.n = n
        self.data_dir = None

    def parse_multifasta(self, file_path):
        """Parse a multifasta file and populate the sequence map."""
        with open(file_path, "r") as file:
            sequence_id = None
            sequence = []
            for line in file:
                line = line.strip()
                if line.startswith(">"):  # New sequence header
                    if sequence_id:
                        self.sequence_map[sequence_id] = ''.join(sequence)
                    sequence_id = line[1:]  # Extract ID
                    sequence = []
                else:
                    sequence.append(line)
            if sequence_id:
                self.sequence_map[sequence_id] = ''.join(sequence)  # Add last sequence

    def create_ngram_database(self):
        """Create an n-gram database from the sequences."""
        for seq_id, sequence in self.sequence_map.items():
            ngrams = [sequence[i:i+self.n] for i in range(len(sequence) - self.n + 1)]
            for ngram in ngrams:
                self.ngram_database[ngram].append(seq_id)

    def save_to_directory(self, directory_path):
        """Save the n-gram database and multifasta file to the specified directory."""
        os.makedirs(directory_path, exist_ok=True)
        self.data_dir = directory_path
       
        # Save the n-gram database
        with open(os.path.join(directory_path, 'ngram_database.pkl'), 'wb') as file:
            pickle.dump(self.ngram_database, file)
       
        # Save the multifasta file
        with open(os.path.join(directory_path, 'multifasta.fasta'), 'w') as file:
            for seq_id, sequence in self.sequence_map.items():
                file.write(f">{seq_id}\n{sequence}\n")

    def load_from_directory(self, directory_path):
        """Load the n-gram database and multifasta file from the specified directory."""
        self.data_dir = directory_path

        # Load the n-gram database
        with open(os.path.join(directory_path, 'ngram_database.pkl'), 'rb') as file:
            self.ngram_database = pickle.load(file)
       
        # Load the multifasta file
        self.sequence_map = {}
        with open(os.path.join(directory_path, 'multifasta.fasta'), 'r') as file:
            sequence_id = None
            sequence = []
            for line in file:
                line = line.strip()
                if line.startswith(">"):  # New sequence header
                    if sequence_id:
                        self.sequence_map[sequence_id] = ''.join(sequence)
                    sequence_id = line[1:]  # Extract ID
                    sequence = []
                else:
                    sequence.append(line)
            if sequence_id:
                self.sequence_map[sequence_id] = ''.join(sequence)  # Add last sequence

    def search_database(self, query_sequence: str, top_x=5):
        """Search the n-gram database and return the top X hits."""
        query_ngrams = [query_sequence[i:i+self.n] for i in range(len(query_sequence) - self.n + 1)]
        scores = defaultdict(int)

        # Count occurrences of query n-grams in the database
        for ngram in query_ngrams:
            if ngram in self.ngram_database:
                for seq_id in self.ngram_database[ngram]:
                    scores[seq_id] += 1

        # Rank sequences by the score
        ranked_sequences = sorted(scores.items(), key=lambda x: x[1], reverse=True)

        # Get top X hits
        top_hits = [(seq_id, scores[seq_id], self.get_sequence(seq_id)) for seq_id, score in islice(ranked_sequences, top_x)]

        return top_hits

    def get_sequence(self, seq_id):
        """Retrieve a sequence by its ID."""
        return self.sequence_map.get(seq_id, None)
   
class PDBAligner:
    def __init__(self, pdb_file1: str, pdb_file2: str):
        self.pdb_file1 = pdb_file1
        self.pdb_file2 = pdb_file2
        self.structure1 = None
        self.structure2 = None
        self.sequence1 = ""
        self.sequence2 = ""
        self.aligned_structure1 = None
        self.aligned_structure2 = None

    def load_structures(self):
        """Load the PDB structures using BioPython's PDBParser."""
        parser = PDBParser(QUIET=True)
        self.structure1 = parser.get_structure("Structure1", self.pdb_file1)
        self.structure2 = parser.get_structure("Structure2", self.pdb_file2)

    def extract_sequences(self):
        """Extract amino acid sequences from the PDB structures."""
        if self.structure1 is None or self.structure2 is None:
            self.load_structures()

        # Use CaPPBuilder to build peptides from the structures
        ppb = CaPPBuilder()

        # Extract sequences from both structures
        sequence1 = []
        sequence2 = []

        # For each structure, build the peptides and get their sequences
        for pp in ppb.build_peptides(self.structure1):
            sequence1.append(pp.get_sequence())

        for pp in ppb.build_peptides(self.structure2):
            sequence2.append(pp.get_sequence())

        # Join sequences from all chains for both structures
        self.sequence1 = "".join(str(seq) for seq in sequence1)
        self.sequence2 = "".join(str(seq) for seq in sequence2)

        return self.sequence1, self.sequence2

    def align_structures(self):
        """Align two PDB structures based on their alpha carbon atoms (CA)."""
        if self.structure1 is None or self.structure2 is None:
            self.load_structures()

        # Extract CA atoms from both structures
        atoms1 = [atom for atom in self.structure1.get_atoms() if atom.get_id() == "CA"]
        atoms2 = [atom for atom in self.structure2.get_atoms() if atom.get_id() == "CA"]

        # Align structures
        super_imposer = Superimposer()
        super_imposer.set_atoms(atoms1[:min(len(atoms1), len(atoms2))],
                                atoms2[:min(len(atoms1), len(atoms2))])
        super_imposer.apply(self.structure2.get_atoms())  # Apply transformation to structure2

        # Store aligned structures
        self.aligned_structure1 = self.structure1
        self.aligned_structure2 = self.structure2

    def generate_alignment_image(self, output_image_path: str = "aligned_structures.png"):
        """Generate an image of the aligned structures using PyMOL."""
        if self.aligned_structure1 is None or self.aligned_structure2 is None:
            self.align_structures()

        # Start a PyMOL session and load the aligned structures
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(self.pdb_file1, "structure1")
            pymol.cmd.load(self.pdb_file2, "structure2")

            # Align the two structures within PyMOL
            pymol.cmd.align("structure2", "structure1")

            # Visualize the structures
            pymol.cmd.hide("everything", "all")
            pymol.cmd.show("cartoon", "structure1")
            pymol.cmd.show("cartoon", "structure2")
            pymol.cmd.color("cyan", "structure1")
            pymol.cmd.color("magenta", "structure2")

            # Reorient the structures in the scene
            pymol.cmd.orient()  # Reorients the scene to fit both structures nicely within view
            pymol.cmd.zoom()    # Ensures the structures are zoomed and centered in the frame
            # Optionally, you can also use pymol.cmd.center("structure1", "structure2") if specific centering is needed

            # Save the image to a specified file
            pymol.cmd.png(output_image_path, width=800, height=600, dpi=300)

        return output_image_path


    def align_sequences(self, line_length=60):
        """Align the sequences of the two PDB structures and return the alignment result."""
        if not self.sequence1 or not self.sequence2:
            self.extract_sequences()

        # Perform global sequence alignment
        alignments = pairwise2.align.globalxx(self.sequence1, self.sequence2)
        alignment = alignments[0]

        # Split the aligned sequences into lines of fixed length
        seq1_aligned = alignment.seqA
        seq2_aligned = alignment.seqB

        alignment_str = []
        for i in range(0, len(seq1_aligned), line_length):
            alignment_str.append(seq1_aligned[i:i+line_length] + "\n" + seq2_aligned[i:i+line_length] + "\n")

        # Join the alignment chunks into a single string
        return "\n".join(alignment_str)

    def proflex_hits(self, hits):
        """Proflex hit output formatting"""
        # Create HTML rows for each hit
        html_rows = "".join([
            f'<tr><td>{hit[0]}</td><td>{hit[1]}</td></tr>\n' for hit in hits
        ])
        return html_rows

    def write_html_output(self, output_html_path: str, output_hits: list, image_output_path: str = "aligned_structures.png"):
        """Generate HTML output embedding the alignment image and sequence alignment."""
        # Generate the structure alignment image
        image_path = self.generate_alignment_image(output_image_path=image_output_path)

        # Align sequences and generate the alignment string
        sequence_alignment = self.align_sequences()

        # Prepare the HTML content
        html_content = f"""
        <html>
        <head>
            <title>Structure and Sequence Alignment</title>
            <style>
                body {{
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                    margin: 0;
                    padding: 0;
                    background-color: #f4f4f9;
                    color: #333;
                }}
                .container {{
                    width: 90%;
                    max-width: 1200px;
                    margin: 20px auto;
                    padding: 20px;
                    background: #f1f1f1;
                    border-radius: 8px;
                    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
                }}
                h1 {{
                    color: #007BFF;
                    font-size: 2.5em;
                    margin-bottom: 0.5em;
                }}
                h2 {{
                    color: #0056b3;
                    border-bottom: 2px solid #007BFF;
                    padding-bottom: 0.5em;
                    font-size: 1.5em;
                    margin-top: 1em;
                }}
                .image-container {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 20px;
                    margin-bottom: 20px;
                }}
                .image-container img {{
                    width: 100%;
                    max-width: 600px;
                    height: auto;
                    border-radius: 4px;
                    border: 1px solid #ddd;
                }}
                .sequence-container, .new-output-container {{
                    font-family: 'Courier New', Courier, monospace;
                    white-space: pre-wrap;
                    border: 1px solid #ddd;
                    padding: 15px;
                    border-radius: 4px;
                    background: #fafafa;
                    overflow-wrap: break-word;
                    word-wrap: break-word;
                    max-width: 100%;
                    margin-bottom: 20px;
                }}
                .sequence-container {{
                    border-color: #007BFF;
                }}
                .new-output-container {{
                    border-color: #28a745;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>Proflex Search Output</h1>
                <h2>Amino Acid Sequence Alignment</h2>
                <div class="image-container">
                    <pre>{sequence_alignment}</pre>
                    <img src="{image_path}" alt="Aligned Structures"/>
                </div>
                <div class="new-output-container">
                    <h2>Top ProFlex Hits</h2>
                    <pre>{output_hits}</pre>
                </div>
            </div>
        </body>
        </html>
        """

        # Write the HTML file
        with open(output_html_path, "w") as f:
            f.write(html_content)

        print(f"HTML output written to {output_html_path}")


class ProFlexQuery:
    def __init__(self, database_directory):
        self.database = NGramDatabase()  # Initialize the NGram database
        self.database.load_from_directory(database_directory)  # Load database from the directory

    def query_pdb(self, query_pdb_file):
        """
        Query a PDB file against the ProFlex database and perform structural and sequence alignment with the top hit.

        Parameters:
        query_pdb_file (str): Path to the query PDB file.

        Returns:
        None
        """
        # Step 1: Perform NMA and ProFlex encoding on the query PDB file
        proflex = ProFlex()  # Initialize ProFlex
        rmsf_dict = proflex.nma_analysis_to_dict(query_pdb_file)  # Perform NMA analysis and get RMSF values
        query_protein_id = list(rmsf_dict.keys())[0]  # Get the protein ID from the PDB title
        rmsf_values = rmsf_dict[query_protein_id]  # Get the RMSF values

        # Encode the RMSF values to a sequence
        encoded_sequence = proflex.encode_sequence(rmsf_values)
        encoded_sequence = "".join(encoded_sequence)
        print(f"Encoded sequence for {query_protein_id}: {encoded_sequence}")

        # Step 2: Search the ProFlex database for the top hit using the encoded sequence
        top_hits = self.database.search_database(encoded_sequence, top_x=5)  # Get the top hit
        if not top_hits:
            print("No hits found in the database.")
            return

        top_hit_id, score, proflex_sequence = top_hits[0]  # Get the top hit's sequence ID
        print(f"Top hit found: {top_hit_id}")
       

        # Step 3: Fetch the PDB model of the top hit
        top_hit_pdb_file = os.path.join(self.database.data_dir, 'PDB', f"{top_hit_id}")
        print(top_hit_pdb_file)
        if not os.path.exists(top_hit_pdb_file):
            print(f"PDB file for top hit {top_hit_id} not found.")
            return

        # Step 4: Perform structural alignment
        pdb_aligner = PDBAligner(query_pdb_file, top_hit_pdb_file)
        pdb_aligner.align_structures()  # Align the structures

        # Step 5: Perform sequence alignment
        query_sequence, top_hit_sequence = pdb_aligner.extract_sequences()  # Extract sequences
        alignment = pdb_aligner.align_sequences()  # Perform sequence alignment

        # Get proflex best hits

        output_hits = pdb_aligner.proflex_hits(top_hits)

        # Step 6: Generate output (e.g., HTML report with structural and sequence alignment)
        output_html_path = f"{query_protein_id}_vs_{top_hit_id}_alignment.html"
        pdb_aligner.write_html_output(output_html_path, output_hits)
        print(f"Results saved to {output_html_path}")


class ProFlex_Aligner:
    def __init__(self, matrix_path, gap_penalty=-5):
        # Load substitution matrix
        self.matrix_df = pd.read_csv(matrix_path, index_col=0)
        self.sub_matrix = self.matrix_df.to_dict()
        self.gap_penalty = gap_penalty

    def needleman_wunsch(self, seq1, seq2):
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m + 1, n + 1))
        traceback = np.zeros((m + 1, n + 1), dtype="object")

        
        for i in range(1, m + 1):
            score_matrix[i][0] = i * self.gap_penalty
            traceback[i][0] = 'up'
        for j in range(1, n + 1):
            score_matrix[0][j] = j * self.gap_penalty
            traceback[0][j] = 'left'

        # fill matrix for alignment
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = self.sub_matrix.get(seq1[i-1], {}).get(seq2[j-1], self.gap_penalty)
                diag = score_matrix[i-1][j-1] + match_score
                up = score_matrix[i-1][j] + self.gap_penalty
                left = score_matrix[i][j-1] + self.gap_penalty

                max_score = max(diag, up, left)
                score_matrix[i][j] = max_score

                if max_score == diag:
                    traceback[i][j] = 'diag'
                elif max_score == up:
                    traceback[i][j] = 'up'
                else:
                    traceback[i][j] = 'left'

        align1, align2 = '', ''
        i, j = m, n
        while i > 0 or j > 0:
            if traceback[i][j] == 'diag':
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif traceback[i][j] == 'up':
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            elif traceback[i][j] == 'left':
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1

        alignment_score = score_matrix[m][n]
        return align1, align2, alignment_score

    def print_alignment(self, align1, align2):
        """basic CLI output"""
        line1 = f"target  0 {align1} {len(align1)}"
        line2 = "         "
        for a, b in zip(align1, align2):
            if a == b:
                line2 += '|'
            elif a == '-' or b == '-':
                line2 += ' '
            else:
                line2 += '.'
        line3 = f"query   0 {align2} {len(align2)}"
        print(line1)
        print(line2)
        print(line3)

    def plot_alignment_heatmap(self, align1, align2, output_path="alignment_heatmap.png"):
        "visual output incorporating log odds scores"
        n = len(align1)
        scores = []

        for a, b in zip(align1, align2):
            if a == "-" or b == "-":
                scores.append(0)
            else:
                scores.append(self.sub_matrix.get(a, {}).get(b, self.gap_penalty))

        # Build annotation matrix
        annotations = np.array([list(align1), list(align2)])

        cmap = sns.diverging_palette(240, 10, as_cmap=True)

        plt.figure(figsize=(n/2, 1.5))
        sns.heatmap(
            np.tile(scores, (2,1)),  # duplicate scores for two rows
            cmap=cmap,
            cbar=True,
            xticklabels=list(range(1, n+1)),
            yticklabels=["target", "query"],
            linewidths=0.5,
            linecolor="gray",
            annot=annotations,
            fmt="",
            annot_kws={"size": 10, "color": "black"},
            vmin=self.matrix_df.min().min(),
            vmax=self.matrix_df.max().max()
        )

        plt.title("Alignment Score Heatmap")
        plt.xlabel("Alignment Position")
        plt.ylabel("")
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"Alignment heatmap saved to {output_path}")
    
    def run_alignment(self, seq1, seq2):
        align1, align2, score = self.needleman_wunsch(seq1, seq2)
        self.print_alignment(align1, align2)
        print(f"Alignment score: {score}")
        self.plot_alignment_heatmap(align1, align2, "ProFlex_alignment.png")


def main():
    print("Welcome to the ProFlex Toolkit")

if __name__ == '__main__':
    main()
