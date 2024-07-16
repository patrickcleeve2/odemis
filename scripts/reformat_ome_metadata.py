#/usr/bin/env python3
import glob
import os
import shutil
import sys
import os

import numpy as np
from ome_types import to_xml
from ome_types.model import OME, Image, Pixels, Channel, TiffData, Plane, Pixels_DimensionOrder, UnitsLength
from ome_types.model.simple_types import PixelType
import tifffile
import sys

# install the following packages
# python3.10 -m venv venv
# source venv/bin/activate
# pip install numpy==1.21.5 ome_types tifffile pylibtiff==0.6.1

# run the script
# navigate to directory containing venv:
# source venv/bin/activate
# python reformat_ome_metadata.py /path/to/ome-tiff-files


def add_odemis_path():
    """Add the odemis path to the python path"""
    def parse_config(path) -> dict:
        """Parse the odemis config file and return a dict with the config values"""
        
        with open(path) as f:
            config = f.read()

        config = config.split("\n")
        config = [line.split("=") for line in config]
        config = {line[0]: line[1].replace('"', "") for line in config if len(line) == 2}
        return config

    odemis_path = "/etc/odemis.conf"
    config = parse_config(odemis_path)
    # dev version,  release version + pyro4
    paths = [f"{config['DEVPATH']}/odemis/src",  "/usr/lib/python3/dist-packages/"]
    for path in paths:
        sys.path.append(path)
    return paths

paths = add_odemis_path()

from odemis import model
from odemis.util.dataio import open_acquisition
from odemis.util import fluo

DEBUG = False # set to True to manipulate synthetic data

def reformat_ome_metadata(filename: str, new_fn: str) -> None:
    """Open an odemis image and reformat the metadata to OME-XML format (2016-06)"""

    # open the odemis image    
    image_data = open_acquisition(filename)

    # get the image dimensions
    size_c = len(image_data)

    d0 = image_data[0].getData()
    if d0.ndim == 5:
        size_t = d0.shape[-4]
        size_z = d0.shape[-3]
    else:
        size_t = 1
        size_z = 1

    size_y = d0.shape[-2]
    size_x = d0.shape[-1]

    # get channel metadata
    channel_md = []
    for d in image_data:

        iwl = d.metadata[model.MD_IN_WL]
        xwl = fluo.get_one_center(iwl) * 1e9  # in nm
        owl = d.metadata[model.MD_OUT_WL]

        # Use excitation wavelength in case of multiple bands
        ewl = fluo.get_one_center_em(owl, iwl) * 1e9  # in nm

        channel_md.append(
            {"emission": ewl, "excitation": xwl}
        )


    # store image data as contig numpy array
    data = np.zeros((size_c, size_t, size_z, size_y, size_x), dtype=np.uint16)

    for c in range(size_c):
        for t in range(size_t):
            for z in range(size_z):
                d = image_data[c].getData()
                data[c, :, :, :, :] = d

    print(f"Image Dimensions (CTZYX):  {data.shape}")

    # for debugging synthetic data (z-index gets brighter)
    if DEBUG:
        print(f"WARN: DEBUG is turned on, images will be manipulated.")
        if size_z > 1:
            for cidx in range(size_c):
                for idx in range(size_z):
                    data[cidx, 0, idx, :, :] = data[cidx, 0, idx, :, :] * (idx+1) 
        # for debugging synthetic data (channel 1 gets dimmer)
        if size_c > 1:
            data[1, :, :, :, :] = 0.5 * data[1, :, :, :, :]

    # get the metadata (same for all images, 
    # apart from channel metadata which is handled separately)
    pixel_size = d0.metadata[model.MD_PIXEL_SIZE]
    pos = d0.metadata[model.MD_POS]
    name = d0.metadata[model.MD_DESCRIPTION]
    exp_time = d0.metadata[model.MD_EXP_TIME]

    # handle single-channel images
    if size_z == 1:
        pos = pos[0], pos[1], None
        pixel_size = pixel_size[0], pixel_size[1], None


    # data blocks
    ifd = 0
    channels = []
    tiff_data_blocks = []
    planes = []
    for nc in range(size_c):
        
        # channel metadata
        channels.append(
            Channel(
                id=f"Channel:0:{nc}",
                name = name,
                illumination_type="Epifluorescence",
                acquisition_mode="WideField",
                contrast_method="Fluorescence",
                excitation_wavelength=channel_md[nc]["excitation"],
                emission_wavelength=channel_md[nc]["emission"],
                samples_per_pixel=1,

        ))

        for nz in range(size_z):

            # add tiff data block
            tiff_data_blocks.append(
                TiffData(
                    plane_count=1,
                    ifd=ifd,
                    first_z=nz,
                    first_c=nc,
                    first_t=0
                )
            )

            # add plane
            planes.append(
                Plane(
                    the_c=nc,
                    the_t=0,
                    the_z=nz,
                    exposure_time=exp_time,
                    position_x=pos[0],
                    position_y=pos[1],
                    position_z=pos[2],
                    position_x_unit=UnitsLength.METER,
                    position_y_unit=UnitsLength.METER,
                    position_z_unit=UnitsLength.METER,
                )
            )

            ifd += 1

    # Create OME metadata structure
    ome = OME()
    image = Image(
        id="Image:0",
        name=os.path.basename(filename),
        acquisition_date=d0.metadata[model.MD_ACQ_DATE],
        pixels=Pixels(
            id="Pixels:0",
            dimension_order=Pixels_DimensionOrder.XYZCT,
            type=PixelType.UINT16,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            size_c=size_c,
            size_t=size_t,
            physical_size_x=pixel_size[0],
            physical_size_y=pixel_size[1],
            physical_size_z=pixel_size[2],
            physical_size_x_unit=UnitsLength.METER,
            physical_size_y_unit=UnitsLength.METER,
            physical_size_z_unit=UnitsLength.METER,
            channels=channels,
            tiff_data_blocks=tiff_data_blocks,
            planes=planes,
        )
    )
    ome.images.append(image)

    # Convert OME object to XML string
    ome_xml = to_xml(ome)

    # Save the image with OME-XML metadata
    with tifffile.TiffWriter(new_fn, bigtiff=True) as tif:
        for t in range(size_t):
            for c in range(size_c):
                for z in range(size_z):
                    tif.write(data[c, t, z], contiguous=True, metadata={'axes': 'YX'})
        tif.overwrite_description(ome_xml)

    print(f"Image saved to {new_fn} with OME-XML metadata.")

def main():
    # get the path from argv
    if len(sys.argv) < 2:
        print("Usage: reformat_ome_metadata.py <path>")
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

    # create a new directory to store the new metadata images
    new_path = os.path.join(PATH, "new-metadata")
    os.makedirs(new_path, exist_ok=True)
    print(f"Creating new directory for reformatted metadata: {new_path}.")

    print(f"Reformatting metadata from {len(filenames)} OME-TIFF files.")
    for fn in filenames:
        print('-'*80)
        new_basename = os.path.basename(
            fn.replace(".ome.tiff", "-new-metadata.ome.tiff")
        )
        new_fn = os.path.join(PATH, "new-metadata", new_basename)

        # reformat the metadata
        reformat_ome_metadata(fn, new_fn)

    print(f"Reformatted metadata from {len(filenames)} OME-TIFF files.")

if __name__ == "__main__":
    main()

    # remove odemis from sys path
    for path in paths:
        print(f"Remove {path} from sys.path")
        try:
            sys.path.remove(path)
        except:
            print(f"Unable to remove {path} from sys.path")