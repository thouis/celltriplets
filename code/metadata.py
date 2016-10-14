import os
import sqlite3
import random


def generate_image_triplets(root, ignore_compounds=["DMSO", "UNKNOWN", "taxol"], exclude_compounds=[]):
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

    # update ignore compounds with exclusions
    ignore_compounds = list(set(ignore_compounds) | set(exclude_compounds))

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
        # compound.  This compound can be in ignore_compounds, but should not
        # be in exclude_compounds
        while True:
            image_3 = random.choice(plate_to_images[image_to_plate[image_1]])
            comp_3, conc_3 = image_to_compoundconc[image_3]
            if comp_3 != comp_1 and comp_3 not in exclude_compounds:
                break

        print("triplet: {} == {} != {}".format(image_1, image_2, image_3))
        print("   CC: {} == {} != {}".format(image_to_compoundconc[image_1],
                                             image_to_compoundconc[image_2],
                                             image_to_compoundconc[image_3]))
        print("   PL: {} == {} != {}".format(image_to_plate[image_1],
                                             image_to_plate[image_2],
                                             image_to_plate[image_3]))
        print("")

if __name__ == "__main__":
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')

    generate_image_triplets(root)
