import sys
import glob
import imread
import json

import numpy as np

from skimage.filters import gaussian_filter
from skimage.filters.rank import minimum
from skimage.morphology import opening, disk
from skimage.feature import peak_local_max


def generate_nuclei_centers(dna):
    # blur a bit
    dna = gaussian_filter(dna, sigma=5) + np.random.uniform(0, 0.01, dna.shape)
    bg = opening(dna, disk(31))
    fgmask = minimum(dna > 2 * bg, disk(5))
    return peak_local_max(dna, min_distance=15, indices=True, labels=fgmask)

if __name__ == '__main__':
    files = glob.glob('{}/*s?_w1*.tif'.format(sys.argv[1]))
    print(len(files))

    centers = {}

    for idx, fname in enumerate(files):
        print(idx, '/', len(files))
        centers[fname] = generate_nuclei_centers(imread.imread(fname)).tolist()

    json.dump(centers, open('{}_nuclei_centers.json'.format(sys.argv[1]), 'w'))
