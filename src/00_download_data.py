import os
import urllib.request
import ssl
import certifi

# Define paths
DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

# Create SSL context with certifi certificates
ssl_context = ssl.create_default_context(cafile=certifi.where())

# Headers to mimic a browser request
headers = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

# 10x Genomics Human Heart URLs
urls = {
    "scRNA.tar.gz": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v3/heart_1k_v3_filtered_feature_bc_matrix.tar.gz",
    "spatial.tar.gz": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Heart/V1_Human_Heart_filtered_feature_bc_matrix.tar.gz",
    "image.tar.gz": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Heart/V1_Human_Heart_spatial.tar.gz"
}

def download_file(url, output_path):
    print(f"Downloading {url}...")
    try:
        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, context=ssl_context) as response:
            with open(output_path, 'wb') as out_file:
                out_file.write(response.read())
        print(f"Saved to {output_path}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")

if __name__ == "__main__":
    print("Starting Download...")
    for filename, url in urls.items():
        output_path = os.path.join(DATA_DIR, filename)
        if not os.path.exists(output_path):
            download_file(url, output_path)
        else:
            print(f"{filename} exists. Skipping.")
    print("Done!")

