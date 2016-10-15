import json
import os
import glob
import imread
import random
import sqlite3


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

        print("triplet: {} == {} != {}".format(image_1, image_2, image_3))
        print("   CC: {} == {} != {}".format(image_to_compoundconc[image_1],
                                             image_to_compoundconc[image_2],
                                             image_to_compoundconc[image_3]))
        print("   PL: {} == {} != {}".format(image_to_plate[image_1],
                                             image_to_plate[image_2],
                                             image_to_plate[image_3]))
        print("")

        yield ("{}/{}".format(image_to_plate[image_1], image_1),
               "{}/{}".format(image_to_plate[image_2], image_2),
               "{}/{}".format(image_to_plate[image_3], image_3))


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

    # infinite generator, never exits
    for im1, im2, im3 in generate_image_triplets(root):
        # make sure there are actually cells to choose from
        if (not centers[im1]) or (not centers[im2]) or (not centers[im3]):
            continue

        yield (random_cell(im1, centers[im1]),
               random_cell(im2, centers[im2]),
               random_cell(im3, centers[im3]))


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

    import numpy as np
    import matplotlib
    matplotlib.use('WXAgg')
    import pylab
    pylab.gray()

    for im1, im2, im3 in generate_triplets(root):
        pylab.imshow(np.hstack((im1, im2, im3)))
        pylab.show()
