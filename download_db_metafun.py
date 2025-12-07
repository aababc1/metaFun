import os
import glob
import hashlib
import json
import tarfile
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from tqdm import tqdm
import argparse
import subprocess
import colorama
import shutil
import sys
MAX_RETRIES = 10
TIMEOUT = 3600  # 1 hour in seconds
CHUNK_SIZE = 8192 * 1024  # 8MB chunks

def download_hashes_file(base_url, output_dir):
    """Download the hashes.json file."""
    url = base_url + "hashes.json"
    file_path = os.path.join(output_dir, "hashes.json")

    response = requests.get(url, timeout=TIMEOUT)
    response.raise_for_status()
    with open(file_path, 'wb') as f:
        f.write(response.content)
    
    with open(file_path, 'r') as f:
        hash_dict = json.load(f)
    
    return hash_dict

BLUE = '\033[94m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BOLD = '\033[1m'
RESET = '\033[0m'


def verify_hash(file_path, expected_hash):
    """Verify the SHA-256 hash of the given file against the expected hash."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    print(f"Verifying hash for {os.path.basename(file_path)}...")
    file_hash = hashlib.sha256()
    
    with open(file_path, "rb") as f, tqdm(unit='B', unit_scale=True, desc=f"Verifying {os.path.basename(file_path)}") as pbar:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            file_hash.update(chunk)
            pbar.update(len(chunk))
    
    if file_hash.hexdigest() != expected_hash:
        raise Exception(f"Hash mismatch for {file_path}")
    
    print(f"Hash verified for {os.path.basename(file_path)}")
    return file_path

def download_file(url, output_dir, expected_hash, retries=0):
    filename = url.split("/")[-1]
    file_path = os.path.join(output_dir, filename)
    temp_file_path = file_path + '.temp'
    
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    
    if os.path.exists(temp_file_path):
        resume_byte_pos = os.path.getsize(temp_file_path)
        headers = {'Range': f'bytes={resume_byte_pos}-'}
        mode = 'ab'  # Append mode
        print(f"Resuming download from byte position {resume_byte_pos}")
    else:
        resume_byte_pos = 0
        headers = {}
        mode = 'wb'  # Write mode
        print("Starting new download")
    
    try:
        with session.get(url, headers=headers, stream=True, timeout=TIMEOUT) as response:
            response.raise_for_status()
            
            if 'Content-Encoding' in response.headers:
                del response.headers['Content-Encoding']
            
            total_size = int(response.headers.get('content-length', 0)) + resume_byte_pos
            
            if resume_byte_pos == total_size:
                print(f"{filename} is already fully downloaded.")
                os.rename(temp_file_path, file_path)
                return verify_hash(file_path, expected_hash)
            
            with open(temp_file_path, mode) as f, tqdm(
                total=total_size,
                initial=resume_byte_pos,
                unit='B',
                unit_scale=True,
                desc=filename
            ) as progress_bar:
                for chunk in response.raw.stream(CHUNK_SIZE, decode_content=False):
                    if chunk:
                        f.write(chunk)
                        progress_bar.update(len(chunk))
        
        os.rename(temp_file_path, file_path)
        return verify_hash(file_path, expected_hash)
    except Exception as e:
        print(f"Error downloading {filename}: {str(e)}. Retrying...")
        if retries < MAX_RETRIES:
            return download_file(url, output_dir, expected_hash, retries + 1)
        else:
            print(f"Failed to download {filename} after {MAX_RETRIES} attempts.")
            raise


def is_sif_file(filename):
    """Check if the file is a SIF image"""
    return filename.endswith('.sif')

def is_tar_gz(filename):
    return filename.endswith('.tar.gz') or filename.endswith('.tgz')


def get_extraction_record_path(output_dir):
    """file path for extraction record"""
    return os.path.join(output_dir, "extracted_files.json")

# def load_extraction_record(output_dir):
#     """load list of extracted files"""
#     record_path = get_extraction_record_path(output_dir)
#     if os.path.exists(record_path):
#         try:
#             with open(record_path, 'r') as f:
#                 return json.load(f)
#         except:
#             return {}
#     return {}
def load_extraction_record(output_dir):
    """Load the list of extracted files"""
    record_file = os.path.join(output_dir, "extracted_files.json")
    try:
        if os.path.exists(record_file):
            with open(record_file, 'r') as f:
                return json.load(f)
        return {}
    except Exception as e:
        print(f"Error loading extraction record: {e}")
        return {}

def save_extraction_record(output_dir, extracted_files):
    """save list of extracted files"""
    record_path = get_extraction_record_path(output_dir)
    with open(record_path, 'w') as f:
        json.dump(extracted_files, f)


def extract_tar_gz(file_path, output_dir):
    if not is_tar_gz(file_path):
        print(f"Skipping extraction for {os.path.basename(file_path)} as it's not a tar.gz file.")
        return False
        
    print(f"Extracting {os.path.basename(file_path)}...")
    success = False
        
    try:
        # Try extraction with pigz
        pigz_process = subprocess.Popen(['pigz', '-dc', file_path], stdout=subprocess.PIPE)
        tar_process = subprocess.Popen(['tar', '-xf', '-', '-C', output_dir], stdin=pigz_process.stdout)
        
        pigz_process.stdout.close()
        tar_process.communicate()
        
        if tar_process.returncode != 0:
            raise subprocess.CalledProcessError(tar_process.returncode, 'tar')
        
        print(f"Extraction complete for {os.path.basename(file_path)}")
        success = True
    except FileNotFoundError:
        print("pigz not found. Falling back to Python's tarfile module.")
        try:
            fallback_extract_tar_gz(file_path, output_dir)
            print(f"Extraction complete using fallback method for {os.path.basename(file_path)}")
            success = True
        except Exception as e:
            print(f"Fallback extraction failed: {e}")
            success = False
    except subprocess.CalledProcessError as e:
        print(f"Error during extraction: {e}")
        print("Falling back to Python's tarfile module.")
        try:
            fallback_extract_tar_gz(file_path, output_dir)
            print(f"Extraction complete using fallback method for {os.path.basename(file_path)}")
            success = True
        except Exception as e:
            print(f"Fallback extraction also failed: {e}")
            success = False
    except Exception as e:
        print(f"Unexpected error during extraction: {e}")
        success = False
    
    return success

def fallback_extract_tar_gz(file_path, output_dir):
    try:
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(path=output_dir)
        print(f"Extraction complete for {os.path.basename(file_path)} using Python's tarfile module")
        return True
    except Exception as e:
        print(f"Python tarfile extraction failed: {e}")
        return False
# def extract_tar_gz(file_path, output_dir):
#     if is_tar_gz(file_path):
#         print(f"Extracting {os.path.basename(file_path)}...")
#         try:
#             pigz_process = subprocess.Popen(['pigz', '-dc', file_path], stdout=subprocess.PIPE)
#             tar_process = subprocess.Popen(['tar', '-xf', '-', '-C', output_dir], stdin=pigz_process.stdout)
            
#             pigz_process.stdout.close()
#             tar_process.communicate()
            
#             if tar_process.returncode != 0:
#                 raise subprocess.CalledProcessError(tar_process.returncode, 'tar')
            
#             print(f"Extraction complete for {os.path.basename(file_path)}")
#         except FileNotFoundError:
#             print("pigz not found. Falling back to Python's tarfile module.")
#             fallback_extract_tar_gz(file_path, output_dir)
#         except subprocess.CalledProcessError as e:
#             print(f"Error during extraction: {e}")
#             print("Falling back to Python's tarfile module.")
#             fallback_extract_tar_gz(file_path, output_dir)
#     else:
#         print(f"Skipping extraction for {os.path.basename(file_path)} as it's not a tar.gz file.")

# def fallback_extract_tar_gz(file_path, output_dir):
#     import tarfile
#     with tarfile.open(file_path, "r:gz") as tar:
#         tar.extractall(path=output_dir)
#     print(f"Extraction complete for {os.path.basename(file_path)} using Python's tarfile module")


#def extract_tar_gz(file_path, output_dir):
#    if is_tar_gz(file_path):
 #       print(f"Extracting {os.path.basename(file_path)}...")
 #       with tarfile.open(file_path, "r:gz") as tar:
 #           tar.extractall(path=output_dir)
 #       print(f"Extraction complete for {os.path.basename(file_path)}")
 #   else:
  #      print(f"Skipping extraction for {os.path.basename(file_path)} as it's not a tar.gz file.")




def merge_kraken2_parts(output_dir, hash_dict):
    """Merge Kraken2 part files into a single archive and extract it"""
    # Identify Kraken2 part files (part_aa, part_ab, etc.)
    part_files = sorted([f for f in os.listdir(output_dir) 
                 if f.startswith("kraken2_GTDBr220_part_")])
    
    if len(part_files) < 4:
        print(f"Cannot merge Kraken2 parts: Found only {len(part_files)} of 4 required parts")
        return False
    
    # Load extracted_files to avoid re-verifying already verified files
    extracted_files = load_extraction_record(output_dir)
    
    # Check if all part files exist
    for part in part_files:
        part_path = os.path.join(output_dir, part)
        if not os.path.exists(part_path):
            print(f"Missing Kraken2 part: {part}")
            return False
    
    # Set merged file path and hash key (case-sensitive)
    merged_filename = os.path.join(output_dir, "kraken2_gtdbr220.tar.gz")
    hash_key = None
    
    # Find exact key in hash_dict (case-insensitive matching)
    for key in hash_dict.keys():
        if key.lower() == "kraken2_gtdbr220.tar.gz":
            hash_key = key
            break
    
    if not hash_key:
        print("Error: Cannot find hash value for merged Kraken2 file in hash_dict")
        return False
    
    try:
        # Start merging if all part files are present
        print(f"Merging {len(part_files)} Kraken2 part files...")
        
        # Merge part files using cat command
        part_paths = [os.path.join(output_dir, part) for part in part_files]
        with open(merged_filename, 'wb') as outfile:
            for part_path in part_paths:
                print(f"Adding {os.path.basename(part_path)}...")
                with open(part_path, 'rb') as infile:
                    # Copy in chunks for memory efficiency
                    shutil.copyfileobj(infile, outfile, length=8*1024*1024)
        
        print("Kraken2 parts merged successfully.")
        print(f"Verifying merged file hash...")
        
        # Hash verification
        merged_hash = hash_dict[hash_key]
        try:
            verify_hash(merged_filename, merged_hash)
            print(f"Hash verified for merged file")
            
            # Extract tar.gz file
            print(f"Extracting merged Kraken2 database...")
            extraction_success = extract_tar_gz(merged_filename, output_dir)
            
            if extraction_success:
                # Save hash to extraction record
                extracted_files[hash_key] = merged_hash
                save_extraction_record(output_dir, extracted_files)
                print(f"Extraction complete and hash saved")
                #decomporess kraken2 gtd r220 database 

                kraken_db_dir = os.path.join(output_dir, "kraken2_GTDBr220")
                print(f"Decompressing Kraken2 database files in {kraken_db_dir}...")

                if os.path.exists(kraken_db_dir):
                    import gzip
                    gz_files = glob.glob(os.path.join(kraken_db_dir, "*.gz"))
                    for gz_file in gz_files:
                        print(f"Decompressing {os.path.basename(gz_file)}...")
                        output_file = gz_file[:-3]  # Remove .gz extension
                        try:
                            with gzip.open(gz_file, 'rb') as f_in:
                                with open(output_file, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                            os.remove(gz_file)  # Remove compressed file after successful decompression
                            print(f"Successfully decompressed {os.path.basename(gz_file)}")
                        except Exception as e:
                            print(f"Error decompressing {gz_file}: {str(e)}")
                else:
                    print(f"Warning: Kraken2 database directory {kraken_db_dir} not found after extraction")
                    

                # Remove merged tar.gz and part files
                os.remove(merged_filename)
                print(f"Removed merged archive")
                
                for part_path in part_paths:
                    os.remove(part_path)
                    print(f"Removed part file: {os.path.basename(part_path)}")
                
                return True
            else:
                print("Failed to extract merged Kraken2 database")
                return False
                
        except Exception as e:
            print(f"Hash verification failed for merged file: {str(e)}")
            if os.path.exists(merged_filename):
                os.remove(merged_filename)
            return False
            
    except Exception as e:
        print(f"Error during Kraken2 merge process: {str(e)}")
        if os.path.exists(merged_filename):
            os.remove(merged_filename)
        return False

def process_file(base_url, output_dir, sif_dir, filename, expected_hash, extracted_files):
    file_url = base_url + filename
    is_sif = is_sif_file(filename)
    target_dir = sif_dir if is_sif else output_dir
    file_path = os.path.join(target_dir, filename)

    # 1. Already verified files (common for all types)
    if filename in extracted_files and extracted_files[filename] == expected_hash:
        print(f"{GREEN if not is_sif else BLUE}Verified (Cached): {filename}{RESET}")
        return file_path  # Return immediately

    # 2. SIF file specific handling
    if is_sif:
        if os.path.exists(file_path):
            try:
                # First-time verification
                verify_hash(file_path, expected_hash)
                extracted_files[filename] = expected_hash
                save_extraction_record(output_dir, extracted_files)
                print(f"{BLUE}SIF Verified: {filename}{RESET}")
                return file_path
            except Exception as e:
                print(f"{RED}SIF Verification Failed: {filename}{RESET}")
                os.remove(file_path)
        
        # New download
        print(f"{YELLOW}Downloading SIF: {filename}{RESET}")
        downloaded_path = download_file(file_url, sif_dir, expected_hash)
        extracted_files[filename] = expected_hash
        save_extraction_record(output_dir, extracted_files)
        return downloaded_path

    # 3. Regular file handling (DB)
    if os.path.exists(file_path):
        try:
            verify_hash(file_path, expected_hash)
            extracted_files[filename] = expected_hash
            save_extraction_record(output_dir, extracted_files)
            print(f"{GREEN}DB Verified: {filename}{RESET}")
            return file_path
        except Exception as e:
            print(f"{RED}DB Verification Failed: {filename}{RESET}")
            os.remove(file_path)
    

    print(f"{YELLOW}Downloading DB: {filename}{RESET}")
    downloaded_path = download_file(file_url, output_dir, expected_hash)
    if is_tar_gz(filename):
        extraction_success = extract_tar_gz(downloaded_path, output_dir)
        if extraction_success:
            extracted_files[filename] = expected_hash
            save_extraction_record(output_dir, extracted_files)
            print(f"Hash saved for {filename}")
            try:
                os.remove(downloaded_path)
                print(f"Original archive {filename} deleted.")
            except Exception as e:
                print(f"Error deleting {filename}: {e}")
        else:
            print(f"Extraction failed for {filename}. Archive not deleted.")
    return downloaded_path


MODULE_FILE_MAP = {
    'RAWREAD_QC': [
        'metafun_v0.1.sif',
        'host_genome.tar.gz'
    ],
    'WMS_TAXONOMY': [
        'metafun_v0.1.sif',
        'kraken2_GTDBr220_part_aa',
        'kraken2_GTDBr220_part_ab', 
        'kraken2_GTDBr220_part_ac',
        'kraken2_GTDBr220_part_ad',
        'sylph_c200_gtdbr220.syldb.tar.gz'
    ],
    'WMS_TAXONOMY_SYLPH': [
        'metafun_v0.1.sif',
        'sylph_c200_gtdbr220.syldb.tar.gz'
    ],
    'WMS_TAXONOMY_KRAKEN': [
        'metafun_v0.1.sif',
        'kraken2_GTDBr220_part_aa',
        'kraken2_GTDBr220_part_ab', 
        'kraken2_GTDBr220_part_ac',
        'kraken2_GTDBr220_part_ad'
    ],
    'INTERACTIVE_TAXONOMY': [
        'interactive_wms_taxonomy_v02.sif'
    ],
    'BIN_ASSESSMENT': [
        'metafun_v0.1.sif',
        'checkM2.tar.gz',
        'gunc.tar.gz',
        'gtdbr220.tar.gz'
    ],
    'COMPARATIVE_ANNOTATION': [
        'metafun_v0.1.sif',
        'interactive_comparative_annotation_v01.sif',
        'eggNOG5.tar.gz',
        'CARD.tar.gz',
        'VFDB.tar.gz',
        'dbCAN.tar.gz',
        'kofam_2023NOV.tar.gz',
        'KEGG_modules.tar.gz'
    ],
    'INTERACTIVE_COMPARATIVE': [
        'interactive_comparative_annotation_v01.sif'
    ],
    'WMS_FUNCTION': [
        'metafun_v0.1.sif',
        'humann3.tar.gz'
    ],
    'GENOME_SELECTOR': [
        'metafun_v0.1.sif'
    ],
    'WMS_STRAIN': [
        'instrain_wms_strain_v03.sif',
        'gtdbr220.tar.gz'
    ],
    'INTERACTIVE_STRAIN': [
        'instrain_wms_strain_v03.sif'
    ],
    'INTERACTIVE_NETWORK': [
        'interactive_wms_network_v04.sif'
    ]
}

def get_file_modules(filename):
    module_order = list(MODULE_FILE_MAP.keys()) 
    return [
        module for module in module_order 
        if filename in MODULE_FILE_MAP[module]
    ]

def check_file_status(filename, expected_hash, extracted_files, sif_dir, output_dir):
    """Correctly determine file status"""
    if is_sif_file(filename):
        file_path = os.path.join(sif_dir, filename)
        

        if filename in extracted_files and extracted_files[filename] == expected_hash:
            return BLUE, '✓'  
            
        if os.path.exists(file_path):
            try:
                verify_hash(file_path, expected_hash)
                extracted_files[filename] = expected_hash
                save_extraction_record(output_dir, extracted_files)
                return BLUE, '✓'
            except FileNotFoundError:
                return RED, '✗'  
            except Exception as e:
                print(f"Hash verification failed for {filename}: {str(e)}")
                return RED, '✗'  
        # File doesn't exist
        return YELLOW, '•'

    file_path = os.path.join(output_dir, filename)

    if filename in extracted_files:
        if extracted_files[filename] == expected_hash:
            return GREEN, '✓'
        return RED, '✗'
    return YELLOW, '•'

def download_and_process(base_url, output_dir, sif_dir, hash_dict, target_module=None, no_kraken=False):
    if no_kraken:
        original_count = len(hash_dict)
        hash_dict = {k: v for k, v in hash_dict.items() if not k.startswith('kraken2_')}
        removed_count = original_count - len(hash_dict)
        print(f"\n[Kraken2 Skip Mode] {removed_count} Kraken2-related files excluded from download")

            
    module_priority = {
        'RAWREAD_QC': [
            'metafun_v0.1.sif',
            'host_genome.tar.gz'
        ],

        'WMS_TAXONOMY': [
            'metafun_v0.1.sif',
            'kraken2_GTDBr220_part_aa',
            'kraken2_GTDBr220_part_ab',
            'kraken2_GTDBr220_part_ac',
            'kraken2_GTDBr220_part_ad',
            'sylph_c200_gtdbr220.syldb.tar.gz'
        ],
        'INTERACTIVE_TAXONOMY': [
            'interactive_wms_taxonomy_v02.sif'
        ],
        'BIN_ASSESSMENT': [
            'metafun_v0.1.sif',
            'checkM2.tar.gz',
            'gunc.tar.gz',
            'gtdbr220.tar.gz'
        ],
        'COMPARATIVE_ANNOTATION': [
            'metafun_v0.1.sif',
            'interactive_comparative_annotation_v01.sif',
            'prokka_latest.sif',
            'eggNOG5.tar.gz',
            'CARD.tar.gz',
            'VFDB.tar.gz',
            'dbCAN.tar.gz',
            'kofam_2023NOV.tar.gz',
            'KEGG_modules.tar.gz'
        ],

        'INTERACTIVE_COMPARATIVE': [
            'interactive_comparative_annotation_v01.sif'
        ],
        'WMS_FUNCTION': [
            'metafun_v0.1.sif',
            'humann3.tar.gz'
        ],
    }

  


# target module designation 
    if target_module:
        required = get_required_files_for_module(target_module, no_kraken)  # Added no_kraken argument
        module_files = required['db'] + required['sif']
        priority_order = module_priority.get(target_module, [])
    
    # Filter only files that actually exist
        valid_priority = [f for f in priority_order if f in hash_dict]
        remaining_files = [f for f in module_files if f not in valid_priority and f in hash_dict]
    
        sorted_files = [(f, hash_dict[f]) for f in valid_priority + remaining_files]
    else:
        # Global priority (for all modules)
        global_priority = [
            'metafun_v0.1.sif',
            'interactive_wms_taxonomy_v02.sif',
            'interactive_comparative_annotation_v01.sif',
            'instrain_wms_strain_v03.sif',
            'interactive_wms_network_v04.sif',
            'host_genome.tar.gz',
            'checkM2.tar.gz',
            'gunc.tar.gz',
            'gtdbr220.tar.gz',
            'sylph_c200_gtdbr220.syldb.tar.gz',
            'humann3.tar.gz',
            'eggNOG5.tar.gz',
            'CARD.tar.gz',
            'VFDB.tar.gz',
            'dbCAN.tar.gz',
            'kofam_2023NOV.tar.gz',
            'KEGG_modules.tar.gz'
        ]
        sorted_files = []
        for f in global_priority:
            if f in hash_dict:
                sorted_files.append((f, hash_dict[f]))
        # Add remaining files
        for f, h in hash_dict.items():
            if f not in global_priority:
                sorted_files.append((f, h))        



    if no_kraken:
        filtered_files = []
        for f in sorted_files:
            if not f[0].startswith("kraken2_"):
                filtered_files.append(f)
            else:
                print(f"Skipping Kraken2 file: {f[0]}")
        sorted_files = filtered_files
  

# load list of extracted files
    extracted_files = load_extraction_record(output_dir)
    sif_files = {}
    kraken_files = {}
    other_files = {}

    

    for filename, hash_value in hash_dict.items():

        if is_sif_file(filename):
            sif_files[filename] = hash_value
        elif filename.startswith("kraken2_"):
            kraken_files[filename] = hash_value
        else:
            other_files[filename] = hash_value
    print("\n--- Checking File Status ---")
    for filename, expected_hash in hash_dict.items():

        color, symbol = check_file_status(filename, expected_hash, extracted_files, sif_dir, output_dir)
        modules = get_file_modules(filename)
        module_info = f" {YELLOW}->{RESET} {', '.join(modules)}" if modules else ""
        
        print(f"{color}{symbol} {filename}{RESET}{module_info}")
    print("=======================\n")

    # 1. sif files download 
    if sif_files:
        print("\n=== Downloading SIF image files ===")
    for filename, expected_hash in sif_files.items():
        try:
            file_path = process_file(base_url, output_dir, sif_dir, filename, expected_hash, extracted_files)
            if file_path:
                print(f"Successfully downloaded {filename} to {file_path}")
            else:
                print(f"File {filename} already verified, skipping download")
        except Exception as e:
            print(f"Failed to download {filename}: {str(e)}")
    if other_files:
        print("\n=== Downloading database files ===")
        for filename, expected_hash in other_files.items():
            try:
                file_path = process_file(base_url, output_dir, sif_dir, filename, expected_hash, extracted_files)
                
                # skip downloaded files 
                if file_path is None:
                    continue
                
                # tar.gz extract
                if file_path and os.path.exists(file_path) and is_tar_gz(filename):
                    print(f"Extracting {filename}...")
                    extraction_success = extract_tar_gz(file_path, output_dir)
                    
                    if extraction_success:
                        # hash save 
                        extracted_files[filename] = expected_hash
                        save_extraction_record(output_dir, extracted_files)
                        print(f"Hash saved for {filename}")
                        
                        # remove tar.gz 
                        try:
                            os.remove(file_path)
                            print(f"Original archive {filename} deleted.")
                        except Exception as e:
                            print(f"Error deleting {filename}: {e}")
                    else:
                        print(f"Extraction failed for {filename}. Archive not deleted.")
                    
            except Exception as e:
                print(f"Failed to process {filename}: {str(e)}")
        print("=== Database downloads complete ===\n")
    
#     if kraken_files and not no_kraken:

# #    if not no_kraken and kraken_files:        
#         print("\n=== Downloading Kraken2 database files (this may take a while) ===")
#         # Kraken2 
#         all_kraken_parts_downloaded = True
#         for filename, expected_hash in kraken_files.items():
#             if filename.lower() == "kraken2_gtdbr220.tar.gz":
#                 continue
                
#             try:
#                 process_file(base_url, output_dir, sif_dir, filename, expected_hash, extracted_files)
#             except Exception as e:
#                 all_kraken_parts_downloaded = False
#                 print(f"Failed to download Kraken2 part {filename}: {str(e)}")

#         print("Merging Kraken2 database parts...")
#         kraken2_parts = [f for f in hash_dict.keys() 
#                         if f.startswith("kraken2_") and 
#                         f != "kraken2_gtdbr220.tar.gz" and
#                         os.path.exists(os.path.join(output_dir, f))]
        
#         if len(kraken2_parts) == 4: 
#             merge_successful = merge_kraken2_parts(output_dir, hash_dict)
#             if merge_successful:
#                 print("Kraken2 database successfully merged and extracted.")
#             else:
#                 print("Failed to merge and extract Kraken2 database.")
#         else:
#             print(f"Cannot merge Kraken2 parts - found only {len(kraken2_parts)} of 4 required parts")
#         print("=== Kraken2 download complete ===\n")
#     #elif kraken_files and no_kraken:
#     #    print("\n=== Skipping Kraken2 database download as requested ===")
  
    if no_kraken:
        print("\n=== Skipping Kraken2 database download as requested ===")
    elif kraken_files:  # Only run when no_kraken is False
        print("\n=== Downloading Kraken2 database files (this may take a while) ===")
        # Kraken2 
        all_kraken_parts_downloaded = True
        for filename, expected_hash in kraken_files.items():
            if filename.lower() == "kraken2_gtdbr220.tar.gz":
                continue
                
            try:
                process_file(base_url, output_dir, sif_dir, filename, expected_hash, extracted_files)
            except Exception as e:
                all_kraken_parts_downloaded = False
                print(f"Failed to download Kraken2 part {filename}: {str(e)}")

        print("Merging Kraken2 database parts...")
        kraken2_parts = [f for f in hash_dict.keys() 
                        if f.startswith("kraken2_") and 
                        f != "kraken2_gtdbr220.tar.gz" and
                        os.path.exists(os.path.join(output_dir, f))]
        
        if len(kraken2_parts) == 4: 
            merge_successful = merge_kraken2_parts(output_dir, hash_dict)
            if merge_successful:
                print("Kraken2 database successfully merged and extracted.")
            else:
                print("Failed to merge and extract Kraken2 database.")
        else:
            print(f"Cannot merge Kraken2 parts - found only {len(kraken2_parts)} of 4 required parts")
        print("=== Kraken2 download complete ===\n")



    return extracted_files 



#get required files for module

#def get_required_files_for_module(module_name, no_kraken=False):  # Added argument
def debug_info(title, message):
    print(f"DEBUG - {title}: {message}")

def get_required_files_for_module(module_name, no_kraken=False):

    """Returns a list of required database and SIF files for a specific module"""
    module_files = {
        "RAWREAD_QC": {
            "db": ["host_genome.tar.gz"],
            "sif": ["metafun_v0.1.sif"]
        },
        "ASSEMBLY_BINNING": {
            "db": [],
            "sif": ["metafun_v0.1.sif"]
        },
        "BIN_ASSESSMENT": {
            "db": ["checkM2.tar.gz", "gunc.tar.gz", "gtdbr220.tar.gz"],
            "sif": ["metafun_v0.1.sif"]
        },
        "WMS_TAXONOMY": {
            "db": ["kraken2_GTDBr220_part_aa", "kraken2_GTDBr220_part_ab", 
                  "kraken2_GTDBr220_part_ac", "kraken2_GTDBr220_part_ad", 
                  "sylph_c200_gtdbr220.syldb.tar.gz"],
            "sif": ["metafun_v0.1.sif"]
        },
        "WMS_TAXONOMY_SYLPH": {
            "db": ["sylph_c200_gtdbr220.syldb.tar.gz"],
            "sif": ["metafun_v0.1.sif"]
        },
        "WMS_TAXONOMY_KRAKEN": {
            "db": ["kraken2_GTDBr220_part_aa", "kraken2_GTDBr220_part_ab", 
                  "kraken2_GTDBr220_part_ac", "kraken2_GTDBr220_part_ad"],
            "sif": ["metafun_v0.1.sif"]
        },
        "WMS_FUNCTION": {
            "db": ["humann3.tar.gz"],
            "sif": ["metafun_v0.1.sif"]
        },
        "COMPARATIVE_ANNOTATION": {
            "db": ["eggNOG5.tar.gz", "CARD.tar.gz", "VFDB.tar.gz", 
                  "dbCAN.tar.gz", "kofam_2023NOV.tar.gz", "KEGG_modules.tar.gz"],
            "sif": ["metafun_v0.1.sif", "interactive_comparative_annotation_v01.sif","prokka_latest.sif"]
        },
        "INTERACTIVE_TAXONOMY": {
            "db": [],
            "sif": ["interactive_wms_taxonomy_v02.sif"]
        },
        "INTERACTIVE_COMPARATIVE": {
            "db": [],
            "sif": ["interactive_comparative_annotation_v01.sif"]
        },
        "GENOME_SELECTOR": {
            "db": [],
            "sif": []
        },
        "WMS_STRAIN": {
            "db": ["gtdbr220.tar.gz"],
            "sif": ["instrain_wms_strain_v03.sif"]
        },
        "INTERACTIVE_STRAIN": {
            "db": [],
            "sif": ["instrain_wms_strain_v03.sif"]
        },
        "INTERACTIVE_NETWORK": {
            "db": [],
            "sif": ["interactive_wms_network_v04.sif"]
        }
    }
    
    for module in module_files:
        if "kraken2_gtdbr220.tar.gz" in module_files[module]["db"]:
            module_files[module]["db"].remove("kraken2_gtdbr220.tar.gz")

    module_name_upper = module_name.upper()
    debug_info("Checking requirements for module", module_name_upper)

    if module_name_upper not in module_files:
        return {"db": [], "sif": []}
    #required = module_files[module_name]
    required = module_files[module_name.upper()]
    debug_info("Required DB files", str(required["db"]))
    debug_info("Required SIF files", str(required["sif"]))
    

    # Remove Kraken2 files if no_kraken is True
    if no_kraken:
        required["db"] = [f for f in required["db"] if not f.startswith("kraken2_")]
    
    return required 


#check files
def check_module_dependencies(module_name, db_dir, sif_dir, no_kraken=False):
    """Check if all required files for a module are available"""
    required = get_required_files_for_module(module_name, no_kraken)
    extraction_record = load_extraction_record(db_dir)
    
    status = {
        "available": True,
        "missing_db": [],
        "missing_sif": []
    }
    
    # Check DB files
    for db_file in required["db"]:
        # Special case for Kraken2 files - check for merged file
        if db_file.startswith("kraken2_"):
            merged_file = "kraken2_gtdbr220.tar.gz"
            if merged_file in extraction_record:
                continue  # Merged Kraken file exists in record
            else:
                status["missing_db"].append(db_file)
                status["available"] = False
                continue
        
        # For tar.gz files - just check extraction record
        if db_file.endswith(".tar.gz"):
            if db_file in extraction_record:
                continue  # File is in extraction record
            else:
                status["missing_db"].append(db_file)
                status["available"] = False
        else:
            # Regular files - check if they exist
            if os.path.exists(os.path.join(db_dir, db_file)):
                continue  # File exists
            else:
                status["missing_db"].append(db_file)
                status["available"] = False
    
    # Check SIF files
    for sif_file in required["sif"]:
        sif_path = os.path.join(sif_dir, sif_file)
        if os.path.exists(sif_path):
            continue  # SIF file exists
        else:
            status["missing_sif"].append(sif_file)
            status["available"] = False
    
    return status

def post_install_verification(output_dir, sif_dir, target_module):
    print("\n=== Post-installation verification started ===")
    verification_passed = True

    # module dependency check
    status = check_module_dependencies(target_module, output_dir, sif_dir)
    if not status['available']:
        print("[CRITICAL] Required files missing:")
        print("\n".join(status['missing_db'] + status['missing_sif']))
        verification_passed = False
    
    # disk space check
    total, used, free = shutil.disk_usage(output_dir)
    if free < 10**9:  # less than 1GB space
        print(f"[warning] disk space is insufficient: {free//10**9}GB")
        verification_passed = False
    
    return verification_passed

def main():
    parser = argparse.ArgumentParser(description="Download metafun databases and SIF images")
    parser.add_argument("--db-dir", help="Directory to store databases", default="./database")
    parser.add_argument("--sif-dir", help="Directory to store SIF images", default="./sif_images")
    parser.add_argument("--module", help="Download files for specific module only") 
    parser.add_argument("--no-kraken", action="store_true", help="Skip downloading Kraken2 database (download later)")
    parser.add_argument("--check-only", action="store_true", help="Only check if files exist, don't download")
    parser.add_argument("--simple-output", action="store_true", help="Simple output with no text, only exit code")
    parser.add_argument("--concise-output", action="store_true", help="Concise output for module dependency check")
    args = parser.parse_args()


    no_kraken = args.no_kraken
    if not no_kraken and not args.check_only:
        print(f"\n{BOLD}WARNING: The Kraken2 database is large (uncompressed: ~400GB){RESET}")
        print(f"Downloading will take a significant amount of time and requires substantial disk space.")
        
        # Get user input (default: n)
        prompt = f"{YELLOW}Do you want to download the Kraken2 database? (y/[n]): {RESET}"
        user_input = input(prompt).strip().lower()

        # Only proceed with Kraken2 download if user enters 'y' or 'yes'
        if user_input in ['y', 'yes']:
            print(f"{GREEN}Proceeding with Kraken2 database download.{RESET}")
            no_kraken = False
        else:
            print(f"{BLUE}Skipping Kraken2 database download.{RESET}")
            no_kraken = True



    if args.check_only and args.module:
        if args.simple_output:
            # Completely quiet mode - return only exit code without any output
            status = check_module_dependencies(args.module, args.db_dir, args.sif_dir, args.no_kraken)
            sys.exit(0 if status["available"] else 1)
        elif args.concise_output:
            status = check_module_dependencies(args.module, args.db_dir, args.sif_dir, args.no_kraken)
            sys.exit(0 if status["available"] else 1)
        else:
            # Concise output mode
            missing = check_module_dependencies(args.module, args.db_dir, args.sif_dir)
            
            color = {
                "RAWREAD_QC": "\033[38;2;255;0;0m",
                "ASSEMBLY_BINNING": "\033[38;2;255;147;0m",
                "BIN_ASSESSMENT": "\033[38;2;0;176;80m",
                "COMPARATIVE_ANNOTATION": "\033[38;2;78;149;217m",
                "WMS_TAXONOMY": "\033[38;2;8;70;250m",
                "WMS_FUNCTION": "\033[38;2;112;48;160m",
                "WMS_STRAIN": "\033[38;2;138;43;226m",
                "INTERACTIVE_COMPARATIVE": "\033[38;2;78;149;217m",
                "INTERACTIVE_TAXONOMY": "\033[38;2;8;70;250m",
                "INTERACTIVE_STRAIN": "\033[38;2;30;144;255m",
                "INTERACTIVE_NETWORK": "\033[38;2;0;206;209m"
            }.get(args.module.upper(), "\033[0m")
            
            # Print status with module color
            print(f"\n{color}--- Module: {args.module} ---\033[0m")

            # DB file status
            for db_file in missing.get("missing_db", []):
                print(f"{YELLOW}• {db_file}\033[0m")

            # SIF file status
            for sif_file in missing.get("missing_sif", []):
                print(f"{YELLOW}• {sif_file}\033[0m")
            
            if missing["missing_db"] or missing["missing_sif"]:
                print(f"\n{RED}Some files are missing for module {args.module}.\033[0m")
                sys.exit(1)
            else:
                print(f"\n{GREEN}All files available for module {args.module}.\033[0m")
                sys.exit(0)

    print(f"\n=== metafun Database Downloader ===")
    print(f"Installing to conda package directories:")
    print(f"Database directory: {os.path.abspath(args.db_dir)}")
    print(f"SIF images directory: {os.path.abspath(args.sif_dir)}")
    
    if not os.access(args.db_dir, os.W_OK):
        raise PermissionError(f"No write permission to database directory: {args.db_dir}")
    if not os.access(args.sif_dir, os.W_OK):
        raise PermissionError(f"No write permission to SIF images directory: {args.sif_dir}")
    
    base_url = "http://www.microbiome.re.kr/home_design/databases/metafun/v1.0/metaFun_db_distribute/"
    #output_dir = "./database"
    
    os.makedirs(args.db_dir, exist_ok=True)
    os.makedirs(args.sif_dir, exist_ok=True)    

    hash_dict = download_hashes_file(base_url, args.db_dir)
    #download_and_process(base_url, args.db_dir, args.sif_dir, hash_dict)
    download_and_process(base_url, args.db_dir, args.sif_dir, hash_dict, args.module, no_kraken) 
    
    print("Download completed successfully!")
    print(f"Databases stored in: {os.path.abspath(args.db_dir)}")
    print(f"SIF images stored in: {os.path.abspath(args.sif_dir)}")


if __name__ == "__main__":
    main()
