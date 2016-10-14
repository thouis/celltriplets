import pylab
import json
import os
import glob
import imread
import random


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

    center_files = glob.glob(os.path.join(root, 'centers', '*json'))
    print("Found {} files of cell centers".format(len(center_files)))

    # combine all centers into one dict
    centers = {}
    for f in center_files:
        [centers.update(json.load(open(f, "r"))) for f in center_files]

    # count the number of nuclei detected
    counts = {k: len(v) for k, v in centers.items()}
    total_centers = sum(counts.values())
    print("{} centers loaded".format(total_centers))

    while True:
        imfile = random.choice(list(centers.keys()))
        im = imread.imread(os.path.join(root, 'images', imfile))
        if counts[imfile]:
            i, j = random.choice(centers[imfile])
            subim = im[max(0, i - 30):, max(0, j - 30):][:60, :60]
            pylab.imshow(subim)
        else:
            pylab.imshow(im)
        pylab.show()
