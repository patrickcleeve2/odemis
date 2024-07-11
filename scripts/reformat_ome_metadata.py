#/usr/bin/env python3
import glob
import os
import shutil
import sys

from libtiff import TIFF
import libtiff.libtiff_ctypes as T

from odemis.util.dataio import open_acquisition
from odemis.dataio import find_fittest_converter
from odemis import model

def reformat_ome_metadata(fn: str, new_fn: str):
    """ Reformat odemis metadata to conform to OME-XML standard. 
    Including Strip ImageJ metadata from an OME-TIFF file
    """

    # copy the file to a new location
    shutil.copy(fn, new_fn)

    # open the image
    data = open_acquisition(fn)
    exporter = find_fittest_converter(new_fn)

    new_data = []
    for i, d in enumerate(data):
        d = d.getData()
        # remove transform
        del d.metadata[model.MD_ROTATION]
        del d.metadata[model.MD_SHEAR]
        # remove extra settings
        del d.metadata[model.MD_EXTRA_SETTINGS]

        # d.metadata["use_um"] = True # -> this is a hack to make the exporter export in micrometers
        # replace the below once the exporter is fixed
        
        # convert the position to micrometers
        # NOTE: this breaks the standard TIFF metadata too... it needs to be fixed in the real exporter
        # QUERY: does the TIFFTAG matter?
        pos = d.metadata[model.MD_POS]
        if len(pos) == 2:
            d.metadata[model.MD_POS] = (pos[0] * 1e6, pos[1] * 1e6)
        elif len(pos) == 3:
            d.metadata[model.MD_POS] = (pos[0] * 1e6, pos[1] * 1e6, pos[2] * 1e6)
        new_data.append(d)

    # save the image
    exporter.export(new_fn, new_data) #only exports a single channel

    # fix x and y position in the metadata (scaled in micrometers)
    old_tfile = TIFF.open(fn, mode="r")
    xpos = old_tfile.GetField(T.TIFFTAG_XPOSITION)
    ypos = old_tfile.GetField(T.TIFFTAG_YPOSITION)
    old_tfile.close()

    # overwrite the metadata inplace
    tfile = TIFF.open(new_fn, mode="r+w")

    # get existing metadata (image description)
    desc = tfile.GetField(T.TIFFTAG_IMAGEDESCRIPTION)

    # strip imagej metadata from the start of the image description
    md_str = desc.decode("utf-8")
    idx = md_str.find("<?xml")

    # set the updated image desc (valid ome-xml)
    tfile.SetField(T.TIFFTAG_IMAGEDESCRIPTION, md_str[idx:].encode("utf-8"))
    
    # set the x and y position for each image # TODO: this only sets for the first page
    tfile.SetField(T.TIFFTAG_XPOSITION, xpos)
    tfile.SetField(T.TIFFTAG_YPOSITION, ypos)

    # save / close file
    tfile.close()

    print(f"reformatted metadata from {fn} to {new_fn}.")


def main():
    # get the path from argv
    if len(sys.argv) < 2:
        print("Usage: strip_metadata.py <path>")
        sys.exit(1)
    PATH = sys.argv[1]

    # check if the path exists
    if not os.path.exists(PATH):
        print(f"Path {PATH} does not exist.")
        sys.exit(1)

    # check if path is directory, if not exit
    if not os.path.isdir(PATH):
        print(f"Path {PATH} is not a directory.")
        sys.exit(1)

    # get all the ome-tiff filenames
    filenames = glob.glob(os.path.join(PATH, "*.ome.tiff"))
    print(f"Found {len(filenames)} OME-TIFF files.")

    # create a new directory to store the stripped metadata images
    new_path = os.path.join(PATH, "stripped-metadata")
    os.makedirs(new_path, exist_ok=True)
    print(f"Creating new directory for stripped metadata: {new_path}.")

    print(f"Stripping metadata from {len(filenames)} OME-TIFF files.")
    for fn in filenames:
        new_basename = os.path.basename(
            fn.replace(".ome.tiff", "-stripped-metadata.ome.tiff")
        )
        new_fn = os.path.join(PATH, "stripped-metadata", new_basename)

        # reformat the metadata
        reformat_ome_metadata(fn, new_fn)

    print(f"Stripped metadata from {len(filenames)} OME-TIFF files.")

if __name__ == "__main__":
    main()
