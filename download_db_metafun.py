import os
import hashlib
import json
import tarfile
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from tqdm import tqdm

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

def is_tar_gz(filename):
    return filename.endswith('.tar.gz') or filename.endswith('.tgz')

def extract_tar_gz(file_path, output_dir):
    if is_tar_gz(file_path):
        print(f"Extracting {os.path.basename(file_path)}...")
        try:
            pigz_process = subprocess.Popen(['pigz', '-dc', file_path], stdout=subprocess.PIPE)
            tar_process = subprocess.Popen(['tar', '-xf', '-', '-C', output_dir], stdin=pigz_process.stdout)
            
            pigz_process.stdout.close()
            tar_process.communicate()
            
            if tar_process.returncode != 0:
                raise subprocess.CalledProcessError(tar_process.returncode, 'tar')
            
            print(f"Extraction complete for {os.path.basename(file_path)}")
        except FileNotFoundError:
            print("pigz not found. Falling back to Python's tarfile module.")
            fallback_extract_tar_gz(file_path, output_dir)
        except subprocess.CalledProcessError as e:
            print(f"Error during extraction: {e}")
            print("Falling back to Python's tarfile module.")
            fallback_extract_tar_gz(file_path, output_dir)
    else:
        print(f"Skipping extraction for {os.path.basename(file_path)} as it's not a tar.gz file.")

def fallback_extract_tar_gz(file_path, output_dir):
    import tarfile
    with tarfile.open(file_path, "r:gz") as tar:
        tar.extractall(path=output_dir)
    print(f"Extraction complete for {os.path.basename(file_path)} using Python's tarfile module")


#def extract_tar_gz(file_path, output_dir):
#    if is_tar_gz(file_path):
 #       print(f"Extracting {os.path.basename(file_path)}...")
 #       with tarfile.open(file_path, "r:gz") as tar:
 #           tar.extractall(path=output_dir)
 #       print(f"Extraction complete for {os.path.basename(file_path)}")
 #   else:
  #      print(f"Skipping extraction for {os.path.basename(file_path)} as it's not a tar.gz file.")

def merge_kraken2_parts(output_dir, hash_dict):
    kraken2_parts = [filename for filename in hash_dict.keys() if filename.startswith("kraken2_") and not is_tar_gz(filename)]
    
    if kraken2_parts:
        print("Merging Kraken2 parts...")
        merged_filename = os.path.join(output_dir, "kraken2_gtdbr220.tar.gz")
        with open(merged_filename, "wb") as outfile:
            for filename in kraken2_parts:
                file_path = os.path.join(output_dir, filename)
                with open(file_path, "rb") as infile:
                    outfile.write(infile.read())
        print("Kraken2 parts merged successfully.")
        
        extract_tar_gz(merged_filename, output_dir)
        os.remove(merged_filename)  # Remove the merged tar.gz file after extraction
        
        for filename in kraken2_parts:
            os.remove(os.path.join(output_dir, filename))

def process_file(base_url, output_dir, filename, expected_hash):
    file_path = os.path.join(output_dir, filename)
    file_url = base_url + filename

    try:
        verify_hash(file_path, expected_hash)
        print(f"{filename} already exists and is verified. Skipping download.")
        return file_path
    except FileNotFoundError:
        print(f"{filename} does not exist. Downloading.")
        return download_file(file_url, output_dir, expected_hash)
    except Exception as e:
        print(f"{filename} hash verification failed: {e}. Re-downloading.")
        if os.path.exists(file_path):
            os.remove(file_path)
        return download_file(file_url, output_dir, expected_hash)

def download_and_process(base_url, output_dir, hash_dict):
    for filename, expected_hash in hash_dict.items():
        try:
            process_file(base_url, output_dir, filename, expected_hash)
        except Exception as e:
            print(f"Failed to process {filename}: {str(e)}")
    
    merge_kraken2_parts(output_dir, hash_dict)
    
    for filename in hash_dict.keys():
        file_path = os.path.join(output_dir, filename)
        if os.path.exists(file_path) and is_tar_gz(file_path) and not filename.startswith("kraken2_"):
            extract_tar_gz(file_path, output_dir)
            os.remove(file_path)  # Remove the tar.gz file after extraction

def main():
    base_url = "http://www.microbiome.re.kr/home_design/databases/metafun/v0.1/metaFun_db_distribute/"
    output_dir = "./database"
    os.makedirs(output_dir, exist_ok=True)
    
    hash_dict = download_hashes_file(base_url, output_dir)
    
    download_and_process(base_url, output_dir, hash_dict)

if __name__ == "__main__":
    main()
