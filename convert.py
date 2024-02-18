from PIL import Image

# source: https://gabrielgambetta.com/computer-graphics-from-scratch/images/crate-texture.jpg
with Image.open('crate-texture.jpg') as im:
    im.save('crate-texture.bmp')