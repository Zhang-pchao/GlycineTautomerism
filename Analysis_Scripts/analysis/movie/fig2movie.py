import imageio
import os
import argparse
import numpy as np
from PIL import Image, ImageDraw, ImageFont

parser = argparse.ArgumentParser(description="Convert all BMP files to an MP4 animation.")
parser.add_argument("--duration", type=float, default=0.1, help="Duration of each frame")
parser.add_argument("--reverse", type=bool, default=False, help="Whether to reverse the animation or not")
parser.add_argument("--loop", type=int, default=1, help="Number of loops, 0 for infinite")
parser.add_argument("--output", type=str, default="out.mp4", help="Output file name")
parser.add_argument("--path", type=str, default=".", help="Path of BMP images")

args = parser.parse_args()

# Create a list of BMP files in the specified path
bmp_files = [f for f in os.listdir(args.path) if f.endswith(".bmp")]

# Sort the BMP files by name
bmp_files.sort()

# Create an Imageio writer object for the MP4 file
writer = imageio.get_writer(args.output, fps=10)

# Initialize the time step
timestep = 0

# Iterate over the BMP files
for file in bmp_files:
    # Load the BMP file as a PIL image
    im = Image.open(os.path.join(args.path, file))
    
    # Create a font object with size 36 points
    font = ImageFont.truetype("/home/pengchao/jupy/ttf/Arial.ttf", 60)
    
    # Create a text label with the current time step
    text = "%4.2f ps" % (timestep * 1e-2)
    
    # Create a text drawing object with black text and white border
    draw = ImageDraw.Draw(im)
    text_width, text_height = draw.textsize(text, font)
    x = im.width - text_width - 10
    y = im.height - text_height - 10
    draw.text((x, y), text, fill="black", font=font, stroke_width=2, stroke_fill="white")
    
    # Convert the PIL image to a NumPy array and append it to the MP4 file
    #im_np = imageio.core.util.array_as_dtype(im, dtype='uint8')
    #writer.append_data(im_np)
    im_np = np.array(im)
    writer.append_data(im_np.astype('uint8'))
    
    # Increment the time step
    timestep += 1
    
    

# Close the Imageio writer object
writer.close()
