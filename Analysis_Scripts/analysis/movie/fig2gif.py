#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

import imageio
import os, sys, glob
import argparse

parser = argparse.ArgumentParser(description = "Convert all bmp files to a gif animation.")

parser.add_argument("--duration", type = float, default = 0.1, help = "Duration of each frame")
parser.add_argument("--reverse", type = bool, default = False, help = "Whether reverse the animation or not")
parser.add_argument("--loop", type = int, default = 1, help = "Loop how many times, 0 for infinite")
parser.add_argument("--output", type = str, default = "out", help = "Whether loop the animation or not")
parser.add_argument("--path", type = str, default = ".", help = "path of .bmp figs")


args = parser.parse_args()

#gif
images = [imageio.imread(filename) for filename in sorted((fn for fn in os.listdir(args.path) if fn.endswith(".bmp")), reverse = args.reverse)]
imageio.mimsave(os.path.join(args.path,args.output) if args.output.endswith(".gif") else ".".join([args.output, "gif"]), images, duration = args.duration, loop =args.loop)

#mp4
#writer = imageio.get_writer('%s.mp4'%args.output, fps=10)
#for file in sorted(os.listdir(args.path)):  
#    if file.endswith(".bmp"):
#        print(file)
#        im = imageio.imread(file)
#        writer.append_data(im)
#writer.close()