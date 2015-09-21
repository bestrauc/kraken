#!/usr/bin/python

import sys
import os
import argparse
import time

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input Basecalls folder to duplicate over time.")
    parser.add_argument("-o", "--output", help="Output directory for BaseCalls folder."
                        " (default: current directory)", default=os.getcwd())
    parser.add_argument("-d", "--delay", help="Delay (s) before writing a new cycle.", 
                        default=120)
    args = parser.parse_args()
    
    return args


def main(argv):
    args = parse_args()
    
    root = os.path.abspath(args.input)

    output_dir = os.path.abspath(args.output) + "/BaseCalls"
    os.mkdir(output_dir)
    
    # descent into the different lane folders
    for lane_dir in os.listdir(root):
        lane_root = os.path.join(root, lane_dir)
        lane_files = os.listdir(lane_root)

        # create the lane folder in the BaseCalls output directory
        lane_copy = os.path.join(output_dir, lane_dir)
        os.mkdir(lane_copy)
        
        # get cycle folders and filter file names
        cycle_dirs = [f for f in lane_files 
                      if os.path.isdir(os.path.join(lane_root, f))]
        filter_files = [f for f in lane_files 
                        if os.path.splitext(f)[1] == ".filter"]

        # sort cyle dirs (by cyle number)..
        cycle_dirs = sorted(cycle_dirs, key = lambda s: int(s.split('.')[0][1:]))
        # .. and filter files (by tile number)
        filter_files = sorted(filter_files, key = lambda s: int(s.split('_')[-1].split('.')[0]))

        #TODO: If we want to write all lanes in parallel, end the loop here,
        # save the paths and encapsulate writing in another loop.

        for f in filter_files:
            os.symlink(os.path.join(lane_root, f), os.path.join(lane_copy, f))

        for d in cycle_dirs:
            os.symlink(os.path.join(lane_root, d), os.path.join(lane_copy, d))
            time.sleep(float(args.delay))
 
        print cycle_dirs
        print filter_files

    
    print argv

if __name__=="__main__":
    main(sys.argv)

