import pandas as pd
import os

def find_extra_fastq_files(metadata_file, fastq_directory):
    """
    Identifies extra fastq.gz files in a directory compared to a metadata file.

    Args:
        metadata_file (str): Path to the metadata CSV file.
        fastq_directory (str): Path to the directory containing fastq.gz files.

    Returns:
        list: A list of filenames that are present in the directory but not expected based on the metadata.
    """

    try:
        # Read the metadata file
        meta_df = pd.read_csv(metadata_file)

        # Create a set of expected fastq.gz filenames
        expected_files = set()
        for new_id in meta_df['new_id']:
            expected_files.add(f"{new_id}_R1.fastq.gz")
            expected_files.add(f"{new_id}_R2.fastq.gz")

        # Get a list of all fastq.gz files in the directory
        actual_files = set(f for f in os.listdir(fastq_directory) if f.endswith(".fastq.gz"))

        # Find the extra files
        extra_files = list(actual_files - expected_files)

        return extra_files

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return []
    except KeyError as e:
        print(f"Error: The 'new_id' column was not found in the metadata file.")
        return []
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return []

# Specify the paths
metadata_file = "~/16s_Arkhangelsk_paper/updated_meta_16s.csv"
fastq_directory = "/mnt/disk1/RUNS/fedorov_de/DAVID/DAVID_BACKUP/draft"

# Find the extra files
extra_files = find_extra_fastq_files(metadata_file, fastq_directory)

# Print the extra files
if extra_files:
    print("The following extra fastq.gz files were found:")
    for filename in extra_files:
        print(filename)
else:
    print("No extra fastq.gz files were found.")
