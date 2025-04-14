import os

def rename_fastq_files(directory, suffix="_B6"):
    """
    Renames FASTQ files in the specified directory by removing technical elements 
    and adding a custom suffix.

    Args:
        directory (str): Path to the directory containing FASTQ files.
        suffix (str): Custom suffix to add to the filenames (default is "_B6").
    """
    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        # Process only files ending with .fastq.gz
        if filename.endswith("001.fastq.gz"):
            # Split the filename into parts
            parts = filename.split("_")
            
            # Extract the sample name (first part) and read number (e.g., R1 or R2)
            sample_name = parts[0]
            read_number = parts[3]  # Assuming "R1" or "R2" is always the 4th part
            
            # Construct the new filename
            new_filename = f"{sample_name}{suffix}_{read_number}.fastq.gz"
            
            # Get full paths for renaming
            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_filename)
            
            # Rename the file
            try:
                os.rename(old_path, new_path)
                print(f"Renamed: {filename} -> {new_filename}")
            except Exception as e:
                print(f"Error renaming file {filename}: {e}")

# Specify the directory containing your FASTQ files
fastq_directory = "/mnt/disk1/RUNS/fedorov_de/DAVID/DAVID_BACKUP/draft"

# Call the function to rename files
rename_fastq_files(fastq_directory)
