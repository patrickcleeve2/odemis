"""
Created on 8 Jan 2025

Copyright © 2025 Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms  of the GNU General Public License version 2 as published by the Free
Software  Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR  PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.
"""

import logging
import sys
from typing import Dict, List, Tuple

import numpy
import yaml

from odemis import model

# tmp: append the path to the 3DCT module until installation solution
sys.path.append("/home/patrick/development/openfibsem/3DCT")
sys.path.append("/home/patrick/development/openfibsem/3DCT/tdct")
from tdct.correlation_v2 import run_correlation
from tdct.util import multi_channel_get_z_guass

def convert_das_to_numpy_stack(das: List[model.DataArray]) -> numpy.ndarray:
    """Convert a list of DataArrays to a numpy stack."""
    arr = []
    for da in das:
        if isinstance(da, model.DataArrayShadow):
            da = da.getData()
        if da.ndim == 5: # model.DataArray is always 5D
            da = da[0, 0, :, :, :]
        arr.append(da)

    return numpy.stack(arr, axis=0)

def get_optimized_z_guass(das: List[model.DataArray], x: int, y: int, z: int, show: bool = False) -> Tuple[float, float, float]:
    """Get the best fitting z-coordinate for the given x, y coordinates. Supports multi-channel images.
    :param das: the data arrays (one per channels)
    :param x: the x-coordinate
    :param y: the y-coordinate
    :param z: the z-coordinate (initial guess)
    :param show: show the plot for debugging
    :return: the z-coordinate and the corresponding correlation value"""
    prev_z = z
    prev_x, prev_y = x, y

    # fm_image  must be 4D np.ndarray with shape (channels, z, y, x)
    fm_image = convert_das_to_numpy_stack(das)

    try:
        # getzGauss can fail, so we need to catch the exception
        zval, z, _ = multi_channel_get_z_guass(image=fm_image, x=x, y=y, show=show)
        logging.info(f"Using Z-Gauss optimisation: {z}, previous z: {prev_z}")

    except RuntimeError as e:
        logging.warning(f"Error in z-gauss optimisation: {e}, using previous z: {prev_z}")
        z = prev_z
        x, y = prev_x, prev_y

    return z


def run_tdct_correlation(fib_coords: numpy.ndarray,
                           fm_coords: numpy.ndarray,
                           poi_coords: numpy.ndarray,
                           fib_image: model.DataArray,
                           fm_image: model.DataArray,
                           path: str) -> Dict[str, Tuple[float, float]]:
    """Run 3DCT Multi-point correlation between FIB and FM images.
    :param fib_coords: the FIB coordinates (x, y, z)
    :param fm_coords: the FM coordinates (x, y, z)
    :param poi_coords: the point of interest coordinates
    :param fib_image: the FIB image
    :param fm_image: the FM image
    :param path: the path to save the results
    :return: the correlation results (including poi, transformations, and errors)
    """

    # get rotation center
    halfmax_dim = int(max(fm_image.shape) * 0.5)
    rotation_center = (halfmax_dim, halfmax_dim, halfmax_dim)

    # get fib pixel size (meters)
    fib_pixel_size = fib_image.metadata[model.MD_PIXEL_SIZE][0]

    # fib image shape minus metadata, fib_pixelsize (microns), fm_image_shape
    image_props = [fib_image.shape, fib_pixel_size * 1e6, fm_image.shape]

    assert fib_image.ndim == 2, "FIB Image must be 2D"
    assert fm_image.ndim == 3, "FM Image must be 3D"
    assert fib_pixel_size is not None, "FIB Pixel Size must be set"
    assert rotation_center is not None, "Rotation Center must be set"
    assert isinstance(rotation_center, tuple), "Rotation Center must be a tuple"
    assert len(rotation_center) == 3, "Rotation Center must have 3 values"

    logging.debug(
        f"Running 3DCT correlation with FIB image shape: {fib_image.shape}, FM image shape: {fm_image.shape}"
    )

    # run correlation
    correlation_results = run_correlation(
        fib_coords=fib_coords,
        fm_coords=fm_coords,
        poi_coords=poi_coords,
        image_props=image_props,
        rotation_center=rotation_center,
        path=path,
        fib_image_filename="FIB-Image",
        fm_image_filename="FM-Image",
    )

    return correlation_results

def get_poi_coordinate(correlation_results: dict) -> Tuple[float, float]:
    """Get the the point of interest coordinate from correlation data
    and convert from micrometers to meters in the microscope image coordinate system.
    :param correlation_results: the correlation results
    :return: the point of interest coordinate in meters
    """
    # get the point of interest coordinate (in microscope coordinates, in metres)
    poi_coord = correlation_results["output"]["poi"][0]["px_um"]
    poi_coord = (poi_coord[0] * 1e-6, poi_coord[1] * 1e-6)
    return poi_coord

def parse_3dct_yaml_file(path: str) -> Tuple[float, float]:
    """Parse the 3DCT yaml file and extract the point of interest (POI)
    in microscope image coordinates (um). Convert the coordinates to metres.
    Note: only the first POI is extracted.
    :param path: Path to the 3DCT yaml file."""
    with open(path, "r") as f:
        data = yaml.safe_load(f)

        pt = get_poi_coordinate(data["correlation"])

    return pt
