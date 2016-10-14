import json
import os
import glob
import imread
import random

from metadata import load_metadata, moa_to_imagelist, image_to_moa, image_to_plate, plate_to_imagelist


def random_cell(imfile, centers, size=60):
    im = imread.imread(os.path.join(root, 'images', imfile))
    i, j = random.choice(centers)
    # make sure full subimage is extracted
    i = min(max(0, i - size // 2), im.shape[0] - size - 1)
    j = min(max(0, j - size // 2), im.shape[1] - size - 1)
    return im[i:(i + size), j:(j + size)]


def generate_triplets(root):
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
        imfile_1 = random.choice(list(centers.keys()))

        moa = image_to_moa[imfile_1]

        # choose image_2 with same moa, but different plate
        while True:
            imfile_2 = random.choice(moa_to_imagelist[moa])
            if image_to_plate[imfile_1] != image_to_plate[imfile_2]:
                break

        # choose image_3 with different moa, but same plate
        while True:
            imfile_3 = random.choice(plate_to_imagelist[image_to_plate[imfile_1]])
            if image_to_moa[imfile_1] == image_to_moa[imfile_3]:
                break

        # make sure there are actually cells to choose from
        if (not centers[imfile_1]) or (not centers[imfile_2]) or (not centers[imfile_3]):
            continue

        yield (random_cell(imfile_1, centers[imfile_1]),
               random_cell(imfile_2, centers[imfile_2]),
               random_cell(imfile_3, centers[imfile_3]))


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

    load_metadata(root)

    for t in generate_triplets(root):
        print(len(t))
