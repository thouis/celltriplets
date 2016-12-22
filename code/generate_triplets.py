import json
import os
import glob
import imread
import random
import sqlite3
import h5py
from collections import defaultdict
import numpy as np


def generate_image_triplets(root, ignore_compounds=["DMSO", "UNKNOWN", "taxol"]):
    conn = sqlite3.connect(os.path.join(root, 'BBBC021_v1.sqlite'))
    c = conn.cursor()

    # get all images, plate names, compounds, and concentrations
    results = c.execute("SELECT Image_Metadata_Plate_DAPI, Image_FileName_DAPI,"
                        "Image_Metadata_Compound, Image_Metadata_Concentration "
                        "from image").fetchall()
    plates, images, compounds, concentrations = zip(*results)
    plate_to_images = {}
    image_to_plate = {}
    for p, i in zip(plates, images):
        plate_to_images.setdefault(p, []).append(i)
        image_to_plate[i] = p

    image_to_compoundconc = {}
    compoundconc_to_images = {}
    for i, cp, cc in zip(images, compounds, concentrations):
        image_to_compoundconc[i] = (cp, cc)
        compoundconc_to_images.setdefault((cp, cc), []).append(i)

    # get all other channels
    results = c.execute("SELECT Image_FileName_DAPI, Image_FileName_Actin, Image_FileName_Tubulin "
                        "from image").fetchall()
    other_images = {d: (a, t) for d, a, t in results}

    while True:
        image_1 = random.choice(images)
        comp_1, conc_1 = image_to_compoundconc[image_1]

        # don't sample from ignored compounds
        if comp_1 in ignore_compounds:
            continue

        # choose a random second image with same compound/conc on a different plate
        while True:
            image_2 = random.choice(compoundconc_to_images[(comp_1, conc_1)])
            if image_to_plate[image_2] != image_to_plate[image_1]:
                break

        # choose a random third image on the same plate, but with a different
        # compound.
        while True:
            image_3 = random.choice(plate_to_images[image_to_plate[image_1]])
            comp_3, conc_3 = image_to_compoundconc[image_3]
            if comp_3 != comp_1 and comp_3 not in ignore_compounds:
                break

        # print("triplet: {} == {} != {}".format(image_1, image_2, image_3))
        # print("   CC: {} == {} != {}".format(image_to_compoundconc[image_1],
        #                                      image_to_compoundconc[image_2],
        #                                      image_to_compoundconc[image_3]))
        # print("   PL: {} == {} != {}".format(image_to_plate[image_1],
        #                                      image_to_plate[image_2],
        #                                      image_to_plate[image_3]))
        # print("")

        yield ((image_to_plate[image_1], image_1) + other_images[image_1],
               (image_to_plate[image_2], image_2) + other_images[image_2],
               (image_to_plate[image_3], image_3) + other_images[image_3])


def cache_image(root, plate, subim, chunksize, _imcache=[None]):
    if _imcache[0] is None:
        _imcache[0] = h5py.File(os.path.join(root, 'image_cache.h5'), 'a')
    cache = _imcache[0]

    grp = cache.require_group(plate)
    if subim not in grp:
        grp.create_dataset(subim,
                           data=imread.imread(os.path.join(root, 'images', plate, subim)),
                           chunks=(chunksize, chunksize))
        cache.flush()
    return grp[subim]


def random_cell(root, iminfo, centers, size=64):
    i, j = random.choice(centers)
    allims = []
    plate = iminfo[0]
    for subim in iminfo[1:]:
        im = cache_image(root, plate, subim, size)
        base_i = min(max(0, i - size // 2), im.shape[0] - size - 1)
        base_j = min(max(0, j - size // 2), im.shape[1] - size - 1)
        subim = im[base_i:(base_i + size), base_j:(base_j + size)]
        allims.append(subim)
    return np.stack(allims, axis=2)


def generate_triplets(root):
    center_files = glob.glob(os.path.join(root, 'centers', '*json'))
    print("Found {} files of cell centers".format(len(center_files)))

    # combine all centers into one dict
    centers = defaultdict(list)
    for f in center_files:
        [centers.update(json.load(open(f, "r"))) for f in center_files]

    # count the number of nuclei detected
    counts = {k: len(v) for k, v in centers.items()}
    total_centers = sum(counts.values())
    print("{} centers loaded".format(total_centers))

    # infinite generator, never exits
    for im1info, im2info, im3info in generate_image_triplets(root):
        im1 = os.path.join(*im1info[:2])
        im2 = os.path.join(*im2info[:2])
        im3 = os.path.join(*im3info[:2])
        # make sure there are actually cells to choose from
        if (not centers[im1]) or (not centers[im2]) or (not centers[im3]):
            continue

        yield (random_cell(root, im1info, centers[im1]),
               random_cell(root, im2info, centers[im2]),
               random_cell(root, im3info, centers[im3]))


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

    import matplotlib
    matplotlib.use('WXAgg')
    import pylab
    pylab.gray()

    for im1, im2, im3 in generate_triplets(root):
        tmp = np.hstack((im1, im2, im3)).astype(np.float32)
        tmp /= tmp.max(axis=(0, 1))
        pylab.imshow(tmp)
        pylab.show()
