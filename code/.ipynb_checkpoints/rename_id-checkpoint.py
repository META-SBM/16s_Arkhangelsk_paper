import pandas as pd
import os
import random

# 1. Read the metadata into a DataFrame
def read_metadata(filepath, id_column):
    """
    Reads metadata from a file and sets the column with fastq.gz filenames as the index.

    Args:
        filepath (str): Path to the metadata file.
        id_column (str): Name of the column with fastq.gz filenames.

    Returns:
        pandas.DataFrame: DataFrame with fastq.gz filenames as the index.
    """
    try:
        df = pd.read_csv(filepath, index_col=id_column)  # Try reading as CSV
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filepath} not found.")
    except ValueError:
        try:
            df = pd.read_table(filepath, index_col=id_column)  # If it fails, try reading as a tab-delimited file
        except FileNotFoundError:
            raise FileNotFoundError(f"File {filepath} not found.")
        except Exception as e:
            raise ValueError(f"Could not read file {filepath}. Check the file format and delimiters.\n{e}")
    return df

# 2. Generate unique IDs
def generate_unique_ids(df, id_prefix="KYH"):
    """
    Generates unique IDs for each sample and adds them to the DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame with samples.
        id_prefix (str): Prefix for unique IDs.

    Returns:
        pandas.DataFrame: DataFrame with an added 'ID_KYH' column.
    """
    n = len(df)
    unique_ids = [f"{id_prefix}{random.randint(100000, 999999)}" for _ in range(n)]
    df['new_id'] = unique_ids
    return df

# 3. Rename fastq files
def rename_fastq_files(df, file_directory=".", extension="_R{read}.fastq.gz"):
    """
    Renames fastq.gz files according to the unique IDs stored in the DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame with 'new_id' column and filenames in the index.
        file_directory (str): Directory where fastq.gz files are located. Defaults to current directory.
        extension (str): File extension template. Defaults to "_R{read}.fastq.gz".
    """
    
    for old_name in df.index:
        for read in [1, 2]:  # For R1 and R2 reads
            new_name = f"{df.loc[old_name, 'new_id']}_R{read}.fastq.gz"
            old_path = os.path.join(file_directory, f'{old_name}_R{read}.fastq.gz')
            new_path = os.path.join(file_directory, new_name)

            try:
                os.rename(old_path, new_path)
                print(f"File {old_name} renamed to {new_name}")
            except FileNotFoundError:
                print(f"File {old_name} not found in directory {file_directory}. Check that the file exists and the path is correct.")
            except OSError as e:
                print(f"Could not rename file {old_name} to {new_name}. You might not have permission or the file already exists. Error: {e}")
    
if __name__ == "__main__":
    # Specify the path to your metadata file
    metadata_file = "~/16s_Arkhangelsk_paper/meta_16s.csv"  # Or "meta_16s.csv" if it's a CSV file

    # Specify the name of the column with fastq.gz filenames
    id_column_name = "df_sample"  # Replace with your actual column name

    # Specify the directory where your fastq.gz files are located (defaults to current directory)
    fastq_directory = "/mnt/disk1/RUNS/fedorov_de/DAVID/DAVID_BACKUP/draft"  

    try:
        # Step 1: Read metadata
        meta_df = read_metadata(metadata_file, id_column_name)

        # Step 2: Generate unique IDs
        meta_df = generate_unique_ids(meta_df)
        # Step 3: Rename fastq files
        rename_fastq_files(meta_df, fastq_directory)

        # (Optional) Save updated table
        meta_df.to_csv("~/16s_Arkhangelsk_paper/updated_meta_16s.csv")
        print("\nUpdated metadata table saved to updated_meta_16s.csv")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")