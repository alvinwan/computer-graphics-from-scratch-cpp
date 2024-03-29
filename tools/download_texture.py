"""
Download and convert texture map used for Raster 12 (Texture)

Before running this file, install prerequisites:

    pip install Pillow requests

Then, run this script from the repository root like the following.

    python tools/download_texture.py
"""

from PIL import Image
import requests

# URL of the texture image
url = "https://gabrielgambetta.com/computer-graphics-from-scratch/images/crate-texture.jpg"

# Local filename to save the image
filename = "./crate-texture.jpg"

# Send a GET request to download the image
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    # Open the file in binary mode for writing
    with open(filename, "wb") as f:
        # Write the image data to the file
        f.write(response.content)
        print(f"Texture downloaded successfully as {filename}")
else:
    print(f"Error downloading texture: {response.status_code}")

# Convert the file from jpg to bmp, a 'raw' format that our C++ code can easily
# parse and use.
outname = './crate-texture.bmp'
with Image.open('crate-texture.jpg') as im:
    im.save(outname)
print(f"Texture converted successfully as {outname}")