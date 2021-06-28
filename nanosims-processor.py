#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import imageio
import numpy as np
import sims
from scipy import ndimage
from skimage import measure
from skimage import transform
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='nanosims processing port')

parser.add_argument('-i', '--input', metavar='PATH',
                    default='GB21_L2_NH4_light_chain1_1.im',
                    dest='input',
                    help='File for input (a .im file)')
parser.add_argument('-f', '--frame', metavar='INT', type=int, default=None,
                    help='Frame index (e.g., "0" would be the first frame" to focus, \
                    default behavior is to loop through all [default: None]')
parser.add_argument('-F', '--filter_method', metavar='STRING', default=None,
                    help='Add a filter step, currently confined to one of "gaussian", "median" [default: None]')
parser.add_argument('-s', '--sigma', metavar='FLOAT', type=float, default=1.0,
                    help='Parameter for filtration (sigma for gaussian, size for median) [default: 1.0]')
parser.add_argument('--rough', action='store_true',
                    help='Only rough plot, then exit')
parser.add_argument('-r', '--roi', metavar='PATH', default='rois.png',
                    help='Path to find manually drawn ROI file (mask) [default: rois.png]')
parser.add_argument('-m', '--res_mult', metavar='PATH', type=float, default=1.65,
                    help='ROI png sometimes gets read in at different dpi than it was saved, \
                    multiplier fixy [default: 1]')
parser.add_argument('-c', '--compare1', metavar='STRING', default='15N 12C',
                    help='Focal element or trolley to be the numerator [default: "15N 12C"]')
parser.add_argument('-C', '--compare2', metavar='STRING', default='14N 12C',
                    help='Comparand element or trolley to be summed with --compare1 for denominator \
                    [default: "14N 12C"]')

args = parser.parse_args()


def apply_filter(img, filt_method, sigma):
    if filt_method == 'gaussian':
        for frame in img.frame.values:
            for specie in img.species.values:
                img.loc[specie, frame,:,:] = ndimage.gaussian_filter(img.loc[specie, frame].data, sigma=sigma)
        return img
    if filt_method == 'median':
        for frame in img.frame.values:
            for specie in img.species.values:
                img.loc[specie, frame, :, :] = ndimage.median_filter(img.loc[specie, frame].data, size=sigma)
    return img

def rough_plot(img):
    for frame in img.frame.values:
        for specie in img.species.values:
            plt.imshow(img.loc[specie, frame], cmap='gray')
            plt.gray()
            plt.axis('off')
            plt.colorbar()
            plt.title(specie)
            plt.show()

def save_plot(img):
    for frame in img.frame.values:
        for specie in img.species.values:
            plt.imshow(img.loc[specie, frame], cmap='gray')
            plt.gray()
            plt.axis('off')
            plt.colorbar()
            plt.title(specie)
            plt.imsave(fname=re.sub(".im$","", args.input) + "_f" + str(frame) + "_" + specie + ".png",
                       arr=aligned_image.loc[specie, frame], cmap='gray')

def get_image_from_raw_rois(roi, im, mult=args.res_mult):

    if roi.shape[0:2] == im.data.shape[2:4]:
        print("The ROI png has the same pixel dimension as the image, great!")
    else:
        print("ROI png dimensions are different from the image, trying to correct...")

        r_h, r_w = [x/2 for x in rois.shape[0:2]]
        im_h, im_w = [x/2 for x in im.data.shape[2:4]]

        #roi = rois[int(r_h - im_h * mult):int(r_h + im_h * mult), int(r_w - im_w * mult):int(r_w + im_w * mult)]

        if r_h == 1366:  # this is half the height of ipad, so hack mode
            not_outside = np.any(roi != [248,248,245], axis=-1)  ## this is color of background
            roi = roi[not_outside,:]

            width_roi = int(np.sqrt(roi.shape[0]))
            width_image = im.data.shape[3]
            res_adj = width_roi / width_image

            roi = roi.reshape((width_roi, width_roi, 3))

            print("just so you know... guessed this is ipad and resolution is off by  "+str(res_adj)+"X so adjusting")

        roi = transform.resize(roi, (width_image, width_image, 3), anti_aliasing=True)*255
    #roi = measure.block_reduce(roi, (round(res_adj), round(res_adj), 1), np.mean)

    return roi

def parse_ROIs(objects, grp_col, c1, c2, annotated_im, stats):
    for obj in range(objects[1]):
        obj_x, obj_y = np.where(objects[0] == (obj + 1))

        counts1 = c1[obj_x, obj_y].sum()
        counts2 = c2[obj_x, obj_y].sum()

        rat_im = counts1/(counts1 + counts2)

        annotated_im[obj_x, obj_y] = rat_im

        stats.append([grp_col, obj, np.array2string(rat_im.values)])

        return annotated_im, stats


image = sims.SIMS(args.input)
aligned_image, shifts = sims.utils.align(image)

# remove 1px border on all sides because Cameca
aligned_image = aligned_image.drop_isel(x=[0, -1], y=[0, -1])

if args.frame:
    aligned_image = aligned_image.loc[:, args.frame,: , :]

if args.filter_method:
    aligned_image = apply_filter(img=aligned_image, filt_method=args.filter_method, sigma=1)

if args.rough:
    rough_plot(aligned_image)
    exit()

try:
    rois = imageio.imread(args.roi, pilmode='RGB')

    # image is in center?
    # hack becasue DPI is different between saved png and imported roi
    rois = get_image_from_raw_rois(roi=rois, im=aligned_image, mult=args.res_mult)

    # set up blank annotated image to be drawn upon
    annotated_image = np.zeros(rois.shape[0:2])

    # 3 groups of rois
    red_rois = rois[:,:,0]-rois[:,:,2]
    red_rois = 1*(red_rois > 200)
    green_rois = rois[:,:,1]-rois[:,:,2]
    green_rois = 1*(green_rois > 200)
    blue_rois = rois[:,:,2]-rois[:,:,1]
    blue_rois = 1*(blue_rois > 200)

    # segment by contiguous colors by category
    red_objects = measure.label(red_rois, return_num=True)
    green_objects = measure.label(green_rois, return_num=True)
    blue_objects = measure.label(blue_rois, return_num=True)

    num_objects = red_objects[1] + green_objects[1] + blue_objects[1]

    stats_table = list()

    for frame in range(aligned_image.data.shape[1]):
        # calculate ratio of desired vs (desired + ref)
        comp1 = aligned_image.loc[args.compare1, frame, :, :]
        comp2 = aligned_image.loc[args.compare2, frame, :, :]

        annotated_image, stats_table = parse_ROIs(objects=red_objects, grp_col='red', c1=comp1, c2=comp2,
                                                  annotated_im=annotated_image, stats=stats_table)
        try:
            annotated_image, stats_table = parse_ROIs(objects=green_objects, grp_col='green', c1=comp1, c2=comp2,
                                                      annotated_im=annotated_image, stats=stats_table)
        except TypeError:
            print("No green, continuing")
        try:
            annotated_image, stats_table = parse_ROIs(objects=blue_objects, grp_col='blue', c1=comp1, c2=comp2,
                                                      annotated_im=annotated_image, stats=stats_table)
        except TypeError:
            print("No blue huh?")

        stats_table = pd.DataFrame(stats_table, columns=["Group", "ROI", "ratio"])

        stats_table.to_csv(re.sub(".im$","", args.input) + "_f0" +
                           "_ratio" + re.sub(" ", "_", args.compare1) +"-x-" + re.sub(" ", "_", args.compare2) +
                           ".tsv", sep="\t", index=False)

        plt.imshow(annotated_image)
        plt.gray()
        plt.axis('off')
        cbar = plt.colorbar()
        cbar.set_label("Fraction (per ROI)")
        plt.title(sims.utils.format_species(args.compare1) + " / (" + sims.utils.format_species(args.compare1)+" + "+
                  sims.utils.format_species(args.compare2) + ")")
        plt.savefig(fname=re.sub(".im$","", args.input) + "_f" + str(frame) +
                          "_ratio" + args.compare1 +"-x-"+ args.compare2 + ".png")
        plt.show()


except FileNotFoundError:
    print('No ROI, saving files for you so you can draw them pls pls ^_^')
    save_plot(aligned_image)

