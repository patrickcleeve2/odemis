# -*- coding: utf-8 -*-
"""
Created on 28 Nov 2013

@author: Kimon Tsitsikas

Copyright © 2012-2013 Kimon Tsitsikas, Delmic

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

from __future__ import division

import numpy
import math
import operator
import scipy.signal
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import transform
import logging
from numpy import unravel_index
from numpy import histogram
from scipy.spatial import cKDTree
from itertools import compress


MAX_STEPS_NUMBER = 100  # How many steps to perform in coordinates matching
SHIFT_THRESHOLD = 0.04  # When to still perform the shift (percentage)
DIFF_NUMBER = 0.95  # Number of values that should be within the allowed difference

def FindCenterCoordinates(subimages):
    """
    For each subimage generated by DivideInNeighborhoods, detects the center 
    of the contained spot. Finally produces a list with the center coordinates 
    corresponding to each subimage.
    subimages (List of model.DataArray): List of 2D arrays containing pixel intensity
    returns (List of tuples): Coordinates of spot centers
    """
    number_of_subimages = len(subimages)
    spot_coordinates = []

    # Pop each subimage from the list
    for i in xrange(number_of_subimages):
        #subimage = subimages[i]
    	# Input might be integer
    	# TODO Dummy, change the way that you handle the array e.g. convolution
    	subimage = subimages[i].astype(numpy.float64)
        subimage_x, subimage_y = subimage.shape

        # See Parthasarathy's paper for details
        xk_onerow = numpy.arange(-(subimage_y - 1) / 2 + 0.5, (subimage_y - 1) / 2, 1)
        (xk_onerow_x,) = xk_onerow.shape
        xk = numpy.tile(xk_onerow, subimage_x - 1)
        xk = xk.reshape((subimage_x - 1, xk_onerow_x))
        yk_onecol = numpy.arange((subimage_x - 1) / 2 - 0.5, -(subimage_x - 1) / 2, -1)
        (yk_onecol_x,) = yk_onecol.shape
        yk_onecol = yk_onecol.reshape((yk_onecol_x, 1))
        yk = numpy.tile(yk_onecol, subimage_y - 1)

        dIdu = subimage[0:subimage_x - 1, 1:subimage_y] - subimage[1:subimage_x, 0:subimage_y - 1]
        dIdv = subimage[0:subimage_x - 1, 0:subimage_y - 1] - subimage[1:subimage_x, 1:subimage_y]

        # Smoothing
        h = numpy.tile(numpy.ones(3) / 9, 3).reshape(3, 3)  # simple 3x3 averaging filter
        dIdu = scipy.signal.convolve2d(dIdu, h, mode='same', fillvalue=0)
        dIdv = scipy.signal.convolve2d(dIdv, h, mode='same', fillvalue=0)

        # Calculate intensity gradient in xy coordinate system
        dIdx = dIdu - dIdv
        dIdy = dIdu + dIdv

        # Assign a,b
        a = -dIdy
        b = dIdx

        # Normalize such that a^2 + b^2 = 1
        I2 = numpy.hypot(a, b)
        s = (I2 != 0)
        a[s] = a[s] / I2[s]
        b[s] = b[s] / I2[s]

        # Solve for c
        c = -a * xk - b * yk

        # Weighting: weight by square of gradient magnitude and inverse distance to gradient intensity centroid.
        dI2 = dIdu * dIdu + dIdv * dIdv
        sdI2 = numpy.sum(dI2[:])
        x0 = numpy.sum(dI2[:] * xk[:]) / sdI2
        y0 = numpy.sum(dI2[:] * yk[:]) / sdI2
        w = dI2 / (0.05 + numpy.sqrt((xk - x0) * (xk - x0) + (yk - y0) * (yk - y0)))

        # Make the edges zero, because of the filter
        w[0, :] = 0
        w[w.shape[0] - 1, :] = 0
        w[:, 0] = 0
        w[:, w.shape[1] - 1] = 0

        # Find radial center
        swa2 = numpy.sum(w[:] * a[:] * a[:])
        swab = numpy.sum(w[:] * a[:] * b[:])
        swb2 = numpy.sum(w[:] * b[:] * b[:])
        swac = numpy.sum(w[:] * a[:] * c[:])
        swbc = numpy.sum(w[:] * b[:] * c[:])
        det = swa2 * swb2 - swab * swab
        xc = (swab * swbc - swb2 * swac) / det
        yc = (swab * swac - swa2 * swbc) / det

        # Output relative to upper left coordinate
        xc = xc + (subimage_y + 1) / 2
        yc = -yc + (subimage_x + 1) / 2
        spot_coordinates.append((xc, yc))

    return spot_coordinates

def DivideInNeighborhoods(data, number_of_spots):
    """
    Given an image that includes N spots, divides it in N subimages with each of them 
    to include one spot. Briefly, it filters the image, finds the N “brightest” spots 
    and crops the region around them generating the subimages. This process is repeated 
    until image division is feasible.
    data (model.DataArray): 2D array containing the intensity of each pixel
    number_of_spots (int,int): The number of CL spots
    returns subimages (List of DataArrays): One subimage per spot
            subimage_coordinates (List of tuples): The coordinates of the center of each 
                                                subimage with respect to the overall image
            subimage_size (int): One dimension because it is square
    """
    subimage_coordinates = []
    subimages = []

    # Filter cosmic ray pixels
    mean_intensity = numpy.mean(data)
    cosmic_ray_thresh = 10 * mean_intensity
    image = numpy.where(data > cosmic_ray_thresh, mean_intensity, data)

    # Determine size of filter window
    filter_window_size = int(image.size / (((number_of_spots[0] * number_of_spots[1]) ** 4))) + 15
    i_max, j_max = unravel_index(image.argmax(), image.shape)
    i_min, j_min = unravel_index(image.argmin(), image.shape)
    max_diff = image[i_max, j_max] - image[i_min, j_min]
    prod_of_spots = numpy.prod(number_of_spots)
    data_max = filters.maximum_filter(image, filter_window_size)
    data_min = filters.minimum_filter(image, filter_window_size)
    

    for i in numpy.arange(3.5, 16.5, 0.1):
        # Determine threshold
        threshold = max_diff / i

        # Filter the parts of the image with variance in intensity greater
        # than the threshold
        maxima = (image == data_max)
        diff = ((data_max - data_min) > threshold)
        maxima[diff == 0] = 0

        labeled, num_objects = ndimage.label(maxima)
        if num_objects > prod_of_spots:
            break

    slices = ndimage.find_objects(labeled)
    
    if len(slices)==0:
        logging.warning("Cannot detect spots.")
        return [],[],0
    
    # Go through these parts and crop the subimages based on the neighborhood_size value
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1) / 2
        y_center = (dy.start + dy.stop - 1) / 2

        subimage_coordinates.append((x_center, y_center))
        # TODO: change +10 and -10 to number relative to spot size
        subimage = image[(dy.start - 10):(dy.stop + 1 + 10), (dx.start - 10):(dx.stop + 1 + 10)]
        if subimage.shape[0]==0 or subimage.shape[1]==0:
            logging.warning("Cannot detect spots.")
            return [],[],0
        subimages.append(subimage)

    # TODO: Handle case where slices is 0 or 1
    # Take care of outliers
    clean_subimages, clean_subimage_coordinates = FilterOutliers(image, subimages, subimage_coordinates)
    subimage_size = subimage.shape[0]

    return clean_subimages, clean_subimage_coordinates, subimage_size

def ReconstructCoordinates(subimage_coordinates, spot_coordinates, subimage_size):
    """
    Given the coordinates of each subimage as also the coordinates of the spot into it, 
    generates the coordinates of the spots with respect to the overall image.
    subimage_coordinates (List of tuples): The coordinates of the 
                                        center of each subimage with 
                                        respect to the overall image
    spot_coordinates (List of tuples): Coordinates of spot centers
    subimage_size(int): One dimension because it is square
    returns (List of tuples): Coordinates of spots in optical image
    """
    optical_coordinates = []
    center_position = (subimage_size / 2) - 1
    for ta, tb in zip(subimage_coordinates, spot_coordinates):
        t = tuple(a + (b - center_position) for a, b in zip(ta, tb))
        optical_coordinates.append(t)

    return optical_coordinates

def FilterOutliers(image, subimages, subimage_coordinates):
    """
    It removes subimages that contain outliers (e.g. cosmic rays).
    image (model.DataArray): 2D array containing the intensity of each pixel
    subimages (List of model.DataArray): List of 2D arrays containing pixel intensity
    returns (List of model.DataArray): List of subimages without the ones containing
                                       outliers
            (List of tuples): The coordinates of the center of each subimage with respect 
                            to the overall image
    """
    number_of_subimages = len(subimages)
    clean_subimages = []
    clean_subimage_coordinates = []
    for i in xrange(number_of_subimages):
        hist, bin_edges = histogram(subimages[i], bins=10)
        # Remove subimage if its histogram implies a cosmic ray
        hist_list = hist.tolist()
        if hist_list.count(0) < 5:
            clean_subimages.append(subimages[i])
            clean_subimage_coordinates.append(subimage_coordinates[i])
            
    # If we removed more than 3 subimages give up and return the initial list
    # This is based on the assumption that each image would contain at maximum
    # 3 cosmic rays.
    if (((len(subimages) - len(clean_subimages)) > 3) or (len(clean_subimages) == 0)):
        clean_subimages = subimages
        clean_subimage_coordinates = subimage_coordinates

    return clean_subimages, clean_subimage_coordinates

def MatchCoordinates(input_coordinates, electron_coordinates, guessing_scale, max_allowed_diff):
    """
    Orders the list of spot coordinates of the grid in the electron image in order to 
    match the corresponding spot coordinates generated by FindCenterCoordinates.
    input_coordinates (List of tuples): Coordinates of spots in optical image
    electron_coordinates (List of tuples): Coordinates of spots in electron image
    guessing_scale (float): Guess scaling for the first transformation
    max_allowed_diff (float): Maximum allowed difference in electron coordinates
    returns (List of tuples): Ordered list of coordinates in electron image with respect 
                                to the order in the electron image
            (List of tuples): List of coordinates in optical image corresponding to the 
                                ordered electron list
    """
    # Remove large outliers
    if len(input_coordinates) > 1:
        optical_coordinates = _FindOutliers(input_coordinates)
    else:
        logging.warning("Cannot find overlay.")
        return [], []

    quality = 0
    guess_x_move, guess_y_move = 0, 0
    guess_scale = guessing_scale
    guess_rotation = 0

    # Informed guess
    guess_coordinates = _TransformCoordinates(optical_coordinates, (guess_x_move, guess_y_move), guess_rotation, (guess_scale, guess_scale))

    # Overlay center
    guess_sub_electron_mean = tuple(map(operator.sub, numpy.mean(guess_coordinates, 0), numpy.mean(electron_coordinates, 0)))
    transformed_coordinates = [tuple(map(operator.sub, guess, guess_sub_electron_mean)) for guess in guess_coordinates]

    max_wrong_points = math.ceil(0.5 * math.sqrt(len(electron_coordinates)))
    for step in xrange(MAX_STEPS_NUMBER):
        #Calculate nearest point
        estimated_coordinates, index1, e_wrong_points, total_shift = _MatchAndCalculate(transformed_coordinates, optical_coordinates, electron_coordinates)

        if estimated_coordinates == []:
            quality = 0
            break

        # Calculate quality
        inv_e_wrong_points = [not i for i in e_wrong_points]
        electron_e_inv_points = [estimated_coordinates[i] for i in list(compress(index1, inv_e_wrong_points))]
        electron_e_points = list(compress(electron_coordinates, inv_e_wrong_points))

        # Calculate distance between the expected and found electron coordinates
        coord_diff = []
        for ta, tb in zip(electron_e_inv_points, electron_e_points):
            tab = tuple(map(operator.sub, ta, tb))
            coord_diff.append(math.hypot(tab[0], tab[1]))

        sort_diff = sorted(coord_diff)
        diff_number_sort = math.floor(DIFF_NUMBER * (len(sort_diff)))
        max_diff = sort_diff[int(diff_number_sort)]

        if max_diff < max_allowed_diff and sum(e_wrong_points) <= max_wrong_points and total_shift <= max_allowed_diff:
            quality = 1
            break

        transformed_coordinates = estimated_coordinates

    if quality == 0:
        logging.warning("Cannot find overlay.")
        return [], []

    # The ordered list gives for each electron coordinate the corresponding optical coordinates
    ordered_coordinates_index = zip(index1, electron_coordinates)
    ordered_coordinates_index.sort()
    ordered_coordinates = []
    for i in xrange(len(ordered_coordinates_index)):
        ordered_coordinates.append(ordered_coordinates_index[i][1])

    # Remove unknown coordinates
    known_ordered_coordinates = list(compress(ordered_coordinates, inv_e_wrong_points))
    known_optical_coordinates = optical_coordinates
    return known_ordered_coordinates, known_optical_coordinates

def _KNNsearch(x_coordinates, y_coordinates):
    """
    Applies K-nearest neighbors search to the lists x_coordinates and y_coordinates.
    x_coordinates (List of tuples): List of coordinates
    y_coordinates (List of tuples): List of coordinates
    returns (List of integers): Contains the index of nearest neighbor in x_coordinates 
                                for the corresponding element in y_coordinates
    """
    points = numpy.array(x_coordinates)
    tree = cKDTree(points)
    distance, index = tree.query(y_coordinates)
    list_index = numpy.array(index).tolist()

    return list_index

def _TransformCoordinates(x_coordinates, translation, rotation, scale):
    """
    Transforms the x_coordinates according to the parameters.
    x_coordinates (List of tuples): List of coordinates
    translation (Tuple of floats): Translation
    rotation (float): Rotation #degrees
    scale (Tuple of floats): Scaling
    returns (List of tuples): Transformed coordinates
    """
    transformed_coordinates = []
    for ta in x_coordinates:
        # translation-scaling-rotation
        translated = tuple(map(operator.add, ta, translation))
        tuple_scale = scale
        scaled = tuple(map(operator.mul, translated, tuple_scale))
        x, y = scaled
        rad_rotation = -rotation * (math.pi / 180)  # rotation in radians, clockwise
        x_rotated = x * math.cos(rad_rotation) - y * math.sin(rad_rotation)
        y_rotated = x * math.sin(rad_rotation) + y * math.cos(rad_rotation)
        rotated = (x_rotated, y_rotated)
        transformed_coordinates.append(rotated)

    return transformed_coordinates

def _MatchAndCalculate(transformed_coordinates, optical_coordinates, electron_coordinates):
    """
    Applies transformation to the optical coordinates in order to match electron coordinates and returns 
    the transformed coordinates. This function must be used recursively until the transformed coordinates
    reach the required accuracy.
    transformed_coordinates (List of tuples): List of transformed coordinates
    optical_coordinates (List of tuples): List of optical coordinates
    electron_coordinates (List of tuples): List of electron coordinates
    returns estimated_coordinates (List of tuples): Estimated optical coordinates
            index1 (List of integers): Indexes of nearest points in optical with respect to electron
            e_wrong_points (List of booleans): Electron coordinates that have no proper match
            total_shift (float): Calculated total shift
    """
    total_shift = 0

    index1 = _KNNsearch(transformed_coordinates, electron_coordinates)
    # Sort optical coordinates based on the _KNNsearch output index
    knn_points1 = [optical_coordinates[i] for i in index1]

    index2 = _KNNsearch(electron_coordinates, transformed_coordinates)
    # Sort electron coordinates based on the _KNNsearch output index
    knn_points2 = [electron_coordinates[i] for i in index2]

    # Sort index1 based on index2 and the opposite
    o_index = [index1[i] for i in index2]
    e_index = [index2[i] for i in index1]

    transformed_range = range(len(transformed_coordinates))
    electron_range = range(len(electron_coordinates))

    # Coordinates that have no proper match (optical and electron)
    o_wrong_points = map(operator.ne, o_index, transformed_range)
    e_wrong_points = map(operator.ne, e_index, electron_range)

    if (all(o_wrong_points) or all(e_wrong_points)):
        logging.warning("Cannot perform matching.")
        return [], [], [], []

    # Calculate the transform parameters for the correct electron_coordinates
    #TODO: Throw exception if inv_e_wrong_points or inv_o_wrong_points has only False elements!!!
    inv_e_wrong_points = [not i for i in e_wrong_points]
    (x_move1, y_move1), (x_scale1, y_scale1), rotation1 = transform.CalculateTransform(list(compress(electron_coordinates, inv_e_wrong_points))
                                 , list(compress(knn_points1, inv_e_wrong_points)))
    x_move1, y_move1 = x_move1, y_move1

    # Calculate the transform parameters for the correct optical_coordinates
    inv_o_wrong_points = [not i for i in o_wrong_points]
    (x_move2, y_move2), (x_scale2, y_scale2), rotation2 = transform.CalculateTransform(list(compress(knn_points2, inv_o_wrong_points))
                                 , list(compress(optical_coordinates, inv_o_wrong_points)))
    x_move2, y_move2 = x_move2, y_move2

    # Average between the two parameters
    #TODO: use numpy.mean()
    avg_x_move = (x_move1 + x_move2) / 2
    avg_y_move = (y_move1 + y_move2) / 2
    avg_x_scale = (x_scale1 + x_scale2) / 2
    avg_y_scale = (x_scale1 + y_scale2) / 2
    avg_rotation = (rotation1 + rotation2) / 2

    # Correct for shift if more than SHIFT_THRESHOLD (percentage) of points are wrong
    # threshold = 2 * SHIFT_THRESHOLD * electron_coordinates.__len__()
    threshold = math.ceil(0.5 * math.sqrt(len(electron_coordinates)))
    # If the number of wrong points is above threshold perform corrections
    if sum(o_wrong_points)>threshold and sum(e_wrong_points)>threshold:
        # Shift
        electron_o_index2 = [electron_coordinates[i] for i in list(compress(index2, o_wrong_points))]
        transformed_o_points = list(compress(transformed_coordinates, o_wrong_points))
        o_wrong_diff = []
        for ta, tb in zip(electron_o_index2, transformed_o_points):
            o_wrong_diff.append(map(operator.sub, ta, tb))

        transformed_e_index1 = [transformed_coordinates[i] for i in list(compress(index1, e_wrong_points))]
        electron_e_points = list(compress(electron_coordinates, e_wrong_points))
        e_wrong_diff = []
        for ta, tb in zip(transformed_e_index1, electron_e_points):
            e_wrong_diff.append(map(operator.sub, ta, tb))

        mean_wrong_diff = numpy.mean(e_wrong_diff,0) - numpy.mean(o_wrong_diff,0)
        avg_x_move = avg_x_move - (0.65 * mean_wrong_diff[0]) / avg_x_scale
        avg_y_move = avg_y_move - (0.65 * mean_wrong_diff[1]) / avg_y_scale
        total_shift = math.hypot((0.65 * mean_wrong_diff[0]) / avg_x_scale, (0.65 * mean_wrong_diff[1]) / avg_x_scale)

        # Angle
        # Calculate angle with respect to its center, therefore move points towards center
        electron_coordinates_vs_center = []
        mean_electron_coordinates = numpy.mean(electron_coordinates, 0)
        for ta in electron_coordinates:
            # translation
            translated = tuple(map(operator.sub, ta, mean_electron_coordinates))
            electron_coordinates_vs_center.append(translated)

        transformed_coordinates_vs_center = []
        for tb in transformed_coordinates:
            # translation
            translated = tuple(map(operator.sub, tb, mean_electron_coordinates))
            transformed_coordinates_vs_center.append(translated)

        # Calculate the angle with its center for every point
        angle_vect_electron = numpy.arctan2([float(i[0]) for i in electron_coordinates_vs_center], [float(i[1]) for i in electron_coordinates_vs_center])
        angle_vect_transformed = numpy.arctan2([float(i[0]) for i in transformed_coordinates_vs_center], [float(i[1]) for i in transformed_coordinates_vs_center])

        # Calculate the angle difference for the wrong electron_coordinates
        angle_vect_transformed_e_index1 = [angle_vect_transformed[i] for i in list(compress(index1, e_wrong_points))]
        angle_diff_electron_wrong = [x - y for x, y in zip(list(compress(angle_vect_electron, e_wrong_points)), angle_vect_transformed_e_index1)]
        angle_diff_electron_wrong[angle_diff_electron_wrong > math.pi] = angle_diff_electron_wrong[angle_diff_electron_wrong > math.pi] - 2 * math.pi
        angle_diff_electron_wrong[angle_diff_electron_wrong < -math.pi] = angle_diff_electron_wrong[angle_diff_electron_wrong < -math.pi] + 2 * math.pi

        # Calculate the angle difference for the wrong transformed_coordinates
        angle_vect_electron_o_index2 = [angle_vect_electron[i] for i in list(compress(index2, o_wrong_points))]
        angle_diff_transformed_wrong = [x - y for x, y in zip(list(compress(angle_vect_transformed, o_wrong_points)), angle_vect_electron_o_index2)]
        angle_diff_transformed_wrong[angle_diff_transformed_wrong > math.pi] = angle_diff_transformed_wrong[angle_diff_transformed_wrong > math.pi] - 2 * math.pi
        angle_diff_transformed_wrong[angle_diff_transformed_wrong < -math.pi] = angle_diff_transformed_wrong[angle_diff_transformed_wrong < -math.pi] + 2 * math.pi

        # Apply correction
        angle_correction = 0.5 * (numpy.mean(angle_diff_electron_wrong, 0) - numpy.mean(angle_diff_transformed_wrong, 0))
        avg_rotation = avg_rotation + 180 / math.pi * angle_correction

    # Perform transformation
    estimated_coordinates = _TransformCoordinates(optical_coordinates, (avg_x_move, avg_y_move), avg_rotation, (avg_x_scale, avg_y_scale))
    new_index1 = _KNNsearch(estimated_coordinates, electron_coordinates)
    new_index2 = _KNNsearch(electron_coordinates, estimated_coordinates)
    new_e_index = [new_index2[i] for i in new_index1]
    new_e_wrong_points = map(operator.ne, new_e_index, electron_range)
    if (all(new_e_wrong_points) or new_index1.count(new_index1[0]) == len(new_index1)):
        logging.warning("Cannot perform matching.")
        return [], [], [], []

    return estimated_coordinates, new_index1, new_e_wrong_points, total_shift

def _FindOutliers(x_coordinates):
    """
    Removes large outliers from the optical coordinates.
    x_coordinates (List of tuples): List of coordinates
    returns (List of tuples): Coordinates without outliers
    """
    # For each point, search for the 2 closest neighbors
    points = numpy.array(x_coordinates)
    tree = cKDTree(points, 2)
    distance, index = tree.query(x_coordinates, 2)
    list_distance = numpy.array(distance)

    # Keep only the second ones because the first ones are the points themselves
    sorted_distance = sorted(list_distance[:, 1])
    outlier_value = 2 * sorted_distance[int(math.ceil(0.5 * len(sorted_distance)))]
    no_outlier_index = list_distance[:, 1] < outlier_value

    return list(compress(x_coordinates, no_outlier_index))
