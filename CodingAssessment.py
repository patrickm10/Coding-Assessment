# Author: Patrick Mejia
import csv
import os
import requests

input_file = "C:/Users/patri/Downloads/sample_files-3/sample_files/annotate/coordinates_to_annotate.txt"
annotation_file = "C:/Users/patri/Downloads/sample_files-3/sample_files/gtf/hg19_annotations.gtf"


# Function to convert the fastq file into a CSV file with Identifier, Sequence, and Length of sequence
# Parameters: FASTQ file and CSV file
# Returns a CSV file with the corresponding FASTQ values and sequence length
def fastq_to_csv(fastq_file, csv_file):
    with open(fastq_file, "r") as f, open(csv_file, "w", newline="") as c:
        writer = csv.writer(c)
        writer.writerow(["Identifier", "Sequence", "Length", "Big Sequence Total"])
        small_seq = 0
        big_seq = 0
        total = 0
        # Parses through every 4 lines since this is the FASTQ format
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = line.strip()
            elif i % 4 == 1:
                sequence = line.strip()
                total += 1
                if len(sequence) < 30:
                    small_seq += 1
                else:
                    big_seq += 1
                writer.writerow([identifier, sequence, len(sequence), big_seq])
        return big_seq, total


# Function to recursively search through FASTQ files in directory
# Parameters: file directory
# Returns percent of sequences in file that are > 30 nucleotides
def find_sequence_percentage(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq"):
                fastq_file = os.path.join(root, file)
                file_name = os.path.basename(fastq_file)
                csv_file = os.path.splitext(fastq_file)[0] + ".csv"
                big_seq, total = fastq_to_csv(fastq_file, csv_file)
                percent = (big_seq / total) * 100
                print(f"File: {file_name}, Sequences longer than 30 nucleotides: {percent}%")


# Function to search for most common sequences in a FASTA file
# Parameters: fasta file
# Returns the 10 most frequent sequences
def find_frequent_sequences(fasta_file):
    sequences = {}
    with open(fasta_file) as f:
        seq = ""
        for line in f:
            if line[0] == ">":
                if seq:
                    # Compares current sequence to sequences in file
                    if seq in sequences:
                        sequences[seq] += 1
                    else:
                        sequences[seq] = 1
                seq = ""
            else:
                seq += line.strip()
        if seq in sequences:
            sequences[seq] += 1
        else:
            sequences[seq] = 1
    freq_seq = sorted(sequences.items(), key=lambda x: x[1], reverse=True)[:10]
    return freq_seq


# Function that parses through the gtf file to get the chromosome name, start and end
# Coordinates and creates a CSV file that improves the search time of the lookup annotation function
# Parameters: gtf file, CSV file
# Returns the CSV file with the corresponding gtf values
def convert_gtf_to_csv(gtf_file, csv_file):
    with open(gtf_file, 'r') as f_in, open(csv_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out)
        writer.writerow(['chromosome', 'start', 'end', 'annotation'])  # write header row
        for line in reader:
            if line[0].startswith('#'):
                continue
            chromosome, start, end, annotation = line[0], int(line[3]), int(line[4]), line[8]
            writer.writerow([chromosome, start, end, annotation])


# Function that takes an input file and a gtf file and looks up the chromosome and coordinate
# If the chromosome exists, it returns the annotated version of the chromosome
# If the chromosome is not found, it returns annotation not found.
# Parameters: chromosome name, starting coordinate, gtf file
# Returns the annotation or a string "Annotation not found"
def lookup_annotation(chromosome, coordinate, gtf_file):
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == chromosome and int(fields[3]) <= coordinate <= int(fields[4]):
                return fields[8]
    return 'Annotation not found'


# Function that reads through each row in an input file and once the lookup function is complete
# it prints the annotation of the chromosome that was found
def find_annotation():
    with open(input_file, 'r') as f:
        for line in f:
            chromosome, coordinate = line.strip().split()
            coordinate = int(coordinate)
            annotation = lookup_annotation(chromosome, coordinate, annotation_file)
            print(f"{chromosome}:{coordinate}: {annotation}")


# Function to calculate the mean target coverage for intervals in a hybrid capture panel
# Parameters: file
# Returns: Dictionary with GC bins in intervals of 10 with their mean target coverages
def mean_target_coverage(file):
    gc_bins = {}
    for i in range(0, 11):
        gc_bins[i * 10] = []

    with open(file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            gc = line[5]
            coverage = line[6]
            try:
                gc = int(float(gc) * 100)
                if 0 <= gc <= 100:
                    bin = gc // 10 * 10
                    coverage = float(coverage)
                    gc_bins[bin].append(coverage)
            # Fixes error when the values cannot be converted
            except ValueError:
                pass
    # Creates dictionary for the final result
    result = {}
    for bin in gc_bins:
        # Checks to find bin values and if so, calculates the mean
        if len(gc_bins[bin]) > 0:
            mean = sum(gc_bins[bin]) / len(gc_bins[bin])
            result[f"{bin} to {bin + 10}"] = mean
    return result


def get_variant_info(ids):
    # Ensembl API endpoint for retrieving variant information
    ensembl = "https://rest.ensembl.org/info/analysis/homo_sapiens?"

    results = []
    for id in ids:
        # Try to reach the Ensembl API
        api = requests.get(ensembl.format(id))
        # Create dataframe for parsed JSON
        data = api.json()

        # Creates datatypes for each extraction value
        location = data["seq_region_name"]
        start = data["start"]
        end = data["end"]

        # Creates a list with the extracted data from Ensembl
        results.append((id, location, start, end))

    return results


# Main code to run file
print("Task 1 - FASTQ")
find_sequence_percentage("C:/Users/patri/Downloads/sample_files-3/sample_files/fastq")
print("---------------------------------------------------------------")
frequent_sequences = find_frequent_sequences("C:/Users/patri/Downloads/sample_files-3/sample_files/fasta/sample.fasta")

print("Task 1 - FASTA")
for seq, count in frequent_sequences:
    print(f"Sequence found {count} times", f" ||  Sequence: {seq}")

print("------------------------------------------------------------------")
print("Task 1 - Lookup Annotation")
file = "C:/Users/patri/Downloads/sample_files-3/sample_files/annotate/coordinates_to_annotate.txt"
print("Uncomment line 183 to see clean output.")
print("I would have used pybedtools library, since it is much more efficient, but PyCharm IDE\n"
      "did not let me download the library in order to use the tools associated.\n")
find_annotation()

print("------------------------------------------------------------------")
print("Task 2")
result = mean_target_coverage("C:/Users/patri/Downloads/sample_files-3/sample_files/Example.hs_intervals.txt")
for bin, mean in result.items():
    print(f"GC bin: {bin}, Mean target coverage: {mean}")

print("------------------------------------------------------------------")
print("Task 3")
print("Attempted task 3, but need more practice using Ensembl API. \n"
      "If you're looking for a software engineer that can learn and adapt very quickly \nand work well in a team, I believe "
      "I can be a very valuable member to Dana Farber.\n"
      "Thanks!", "\U0001f600")
