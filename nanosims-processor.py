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
from matplotlib import colors

parser = argparse.ArgumentParser(description = 'Tool for parsing Cameca NanoSIMS .im files to produce ratios. \
                                 Basic usage is to read a .im file, and produce single-channel PNGs and/or ratios of \
                                 isotope_1 / (isotope_1 + isotope_2). Also can read in a .png as a mask for ROIs.')
group1 = parser.add_argument_group('BASIC OPTIONS', '')
group2 = parser.add_argument_group('PLOTTING OPTIONS', '')
group3 = parser.add_argument_group('ROI OPTIONS', '')

group1.add_argument('-i', '--input', metavar='PATH', required=True,
                    dest='input', help='File for input (a .im file)')
group1.add_argument('-o', '--output', metavar='PATH',
                    default=None, help='Prefix for output [default: --input without the .im ending]')
group1.add_argument('-c', '--compare1', metavar='"STRING"', default='15N 12C',
                    help='Focal element or trolley to be the numerator (use quotes if has spaces) [default: "15N 12C"]')
group1.add_argument('-C', '--compare2', metavar='"STRING"', default='14N 12C',
                    help='Comparand element or trolley to be summed with --compare1 for denominator \
                    [default: "14N 12C"]')
group1.add_argument('-f', '--frame', metavar='INT', type=int, default=None,
                    help='Frame index (e.g., "0" would be the first frame) to focus, \
                    default behavior is to loop through all [default: None]')
group1.add_argument('-F', '--filter_method', metavar='STRING', default=None,
                    help='Add a filter step, currently confined to one of "gaussian", "median" [default: None]')
group1.add_argument('-s', '--sigma', metavar='FLOAT', type=float, default=1.0,
                    help='Parameter for filtration (sigma for gaussian, size for median) [default: 1.0]')
group2.add_argument('--rough', action='store_true',
                    help='Only rough plot, then exit [default: False]')
group2.add_argument('-w', '--whole_image', action='store_true',
                    help='Calculate ratio of specified channels for the full image (all pixels)')
group2.add_argument('-n', '--no_show', action='store_true',
                    help='Add this to not show plots to your screen, only save. Useful when batch processing')
group3.add_argument('-r', '--roi', metavar='PATH', default=None,
                    help='Path to corresponding manually drawn ROI file (mask) [default: None]')
group3.add_argument('--count_stat', metavar='STRING', default='mean',
                    help='Method for aggregating pixel values in each ROI to report in the summary table. \
                    Does not affect ratios. Current options: mean, sum [default: mean]')

args = parser.parse_args()

def apply_filter(img, filt_method, sigma):
    if filt_method == 'gaussian':
        for frame in img.frame.values:
            for specie in img.species.values:
                img.loc[specie, frame,: ,:] = ndimage.gaussian_filter(img.loc[specie, frame].data, sigma=sigma)
        return img
    if filt_method == 'median':
        for frame in img.frame.values:
            for specie in img.species.values:
                img.loc[specie, frame, :, :] = ndimage.median_filter(img.loc[specie, frame].data, size=int(sigma))
    return img

def rough_plot(img):
    for frame in img.frame.values:
        for specie in img.species.values:
            plt.imshow(img.loc[specie, frame], cmap='gray', interpolation='none')
            plt.gray()
            plt.axis('off')
            plt.colorbar()
            plt.title(specie)
            plt.show()

def save_plot(img):
    for frame in img.frame.values:
        for specie in img.species.values:
            plt.cla()
            plt.imshow(img.loc[specie, frame], cmap='gray', interpolation='none')
            plt.gray()
            plt.axis('off')
            plt.colorbar()
            plt.title(specie)
            plt.imsave(fname=args.output + "_f" + str(frame) + "_" + specie + ".png",
                       arr=aligned_image.loc[specie, frame], cmap='gray')

def get_image_from_raw_rois(roi, im):

    if roi.shape[0:2] == im.data.shape[2:4]:
        print("The ROI png has the same pixel dimension as the image, great! [no ROI cropping]")
    elif roi.shape[0] > im.data.shape[2]:
        print("ROI png dimensions are larger than the image, trying to correct...")

        edges = np.concatenate([roi[[0, -1], 0, :], roi[[0, -1], 0, :]]).reshape((4, 3))
        if np.all(edges == edges[0, :]):  # check if edges are all the same RGB (likely border)
            print("Edge pixels are all the same color, assuming this is background so trimming these solid edges")

            not_outside = np.any(roi != edges[0, :], axis=-1)
            roi = roi[not_outside, :]

            width_roi = int(np.sqrt(roi.shape[0]))
            width_image = im.data.shape[3]

            roi = roi.reshape((width_roi, width_roi, 3))

        # check if this fixed it
        if roi.shape[0:2] == im.data.shape[2:4]:
            print("Cropped successfully, no rescaling!")
        else:
            print("FYI... ROI was still too big @ "+str(width_roi)+"px vs. "+str(width_image)+"px so scaling")
            roi = transform.resize(roi, (width_image, width_image, 3), anti_aliasing=True)*255
    else:
        print("Your ROI is smaller? Stopping so you can fix. Image: " + im.data.shape + " vs roi: " + roi.shape)
        exit()

    return roi

def parse_ROIs(objects, grp_col, c1, c2, annotated_im, im, stats):
    for obj in range(objects[1]):
        obj_x, obj_y = np.where(objects[0] == (obj + 1))

        pts = np.where(objects[0] == (obj + 1))
        pts = np.reshape(pts, (2, len(pts[0])))
        pts = [pts[:, x] for x in range(pts.shape[1])]

        counts1 = sum([c1.data[x, y] for x, y in pts])
        counts2 = sum([c2.data[x, y] for x, y in pts])

        if args.count_stat == 'mean':
            all_counts = [x for x in np.array([im.data[:, x, y] for x, y in pts]).sum(axis=0)]
        elif args.count_stat == 'sum':
            all_counts = [x for x in np.array([im.data[:, x, y] for x, y in pts]).mean(axis=0)]

        rat_im = counts1/(counts1 + counts2)
        annotated_im[obj_x, obj_y] = rat_im

        obj_stats = [grp_col, obj]
        obj_stats.extend(all_counts)
        obj_stats.append(rat_im)

        stats.append(obj_stats)
    return annotated_im, stats

### BODY ###
image = sims.SIMS(args.input)
aligned_image, shifts = sims.utils.align(image)

if args.compare1 not in aligned_image.species:
    print("Your --compare1 value is not in the .im file? Did you put an underscore in? Available names:\n" +
          ", ".join(aligned_image.species.values))
    exit()
if args.compare2 not in aligned_image.species:
    print("Your --compare2 value is not in the .im file? Did you put an underscore in? Available names:\n" +
          ", ".join(aligned_image.species.values))
    exit()

# remove 1px border on all sides because Cameca
aligned_image = aligned_image.drop_isel(x=[0, -1], y=[0, -1])

if args.frame is not None:
    print('subsetting to frame ' + str(args.frame))
    aligned_image = aligned_image.loc[:, [args.frame], : , :]
if args.filter_method:
    aligned_image = apply_filter(img=aligned_image, filt_method=args.filter_method, sigma=args.sigma)
if args.output is None:
    args.output = re.sub(".im", "", args.input)
if args.rough:
    rough_plot(aligned_image)
    exit()

try:
    rois = imageio.imread(args.roi, pilmode='RGB')

    # hack becasue DPI is different between saved png and imported roi
    rois = get_image_from_raw_rois(roi=rois, im=aligned_image)

    # set up blank annotated image to be drawn upon
    annotated_image = np.zeros(rois.shape[0:2])

    # 3 groups of rois - red, green, blue
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

    # make a nice annotated image so rois are numbered to correspond to future stats_table
    obj = np.where(red_objects[0] > 0, red_objects[0], 0) + \
          np.where(green_objects[0] > 0, green_objects[0] + red_objects[1], 0) + \
          np.where(blue_objects[0] > 0, blue_objects[0] + red_objects[1] + green_objects[1], 0)

    cmap, norm = colors.from_levels_and_colors(levels=np.unique(obj), extend='max',
                                               colors=np.hstack((np.repeat('black', 1),
                                                                 np.repeat('red', red_objects[1]),
                                                                 np.repeat('green', green_objects[1]),
                                                                 np.repeat('blue', blue_objects[1]))))

    plt.imshow(obj, cmap=cmap, norm=norm, interpolation='none')
    for region in measure.regionprops(measure.label(red_rois)):
        ry, rx = region.centroid
        plt.text(rx, ry, s=region.label, ha='center', va='center', bbox=dict(fc='white', ec='none', pad=2))
    for region in measure.regionprops(measure.label(blue_rois)):
        ry, rx = region.centroid
        plt.text(rx, ry, s=region.label, ha='center', va='center', bbox=dict(fc='white', ec='none', pad=2))
    for region in measure.regionprops(measure.label(blue_rois)):
        ry, rx = region.centroid
        plt.text(rx, ry, s=region.label, ha='center', va='center', bbox=dict(fc='white', ec='none', pad=2))
    plt.savefig(fname=args.output + "_roi_table.png")
    if not args.no_show:
        plt.show()
    plt.cla()

    stats_table = list()

except OSError:
    if not args.whole_image:
        print('No ROI, saving pngs for each channel so you can use for making ROIs (or rerun with -r)')
        save_plot(aligned_image)


for frame in range(aligned_image.data.shape[1]):
    # calculate ratio of desired vs (desired + ref)
    comp1 = aligned_image.loc[args.compare1, frame, :, :]
    comp2 = aligned_image.loc[args.compare2, frame, :, :]

    if args.whole_image:
        ratio_image = comp1/(comp1 + comp2)
        ratio_image = np.nan_to_num(ratio_image)

        plt.cla()
        plt.imshow(ratio_image, interpolation='none', cmap='cividis')
        plt.axis('off')
        plt.colorbar().set_label("Fraction (per ROI)")
        plt.title(sims.utils.format_species(args.compare1) + " / (" + sims.utils.format_species(args.compare1)+" + "+
                  sims.utils.format_species(args.compare2) + ")")
        plt.savefig(fname=args.output + "_whole_f" + str(frame) +
                          "_ratio" + re.sub(" ", "_", args.compare1) +"-x-"+ re.sub(" ", "_", args.compare2) + ".png")
        if not args.no_show:
            plt.show()

    if args.roi:
        try:
            annotated_image, stats_table = parse_ROIs(objects=red_objects, grp_col='red', c1=comp1, c2=comp2,
                                                      im=aligned_image.loc[:, frame, :, :],
                                                      annotated_im=annotated_image, stats=stats_table)
        except TypeError:
            print("No red, continuing")
        try:
            annotated_image, stats_table = parse_ROIs(objects=green_objects, grp_col='green', c1=comp1, c2=comp2,
                                                      im=aligned_image.loc[:, frame, :, :],
                                                      annotated_im=annotated_image, stats=stats_table)
        except TypeError:
            print("No green, continuing")
        try:
            annotated_image, stats_table = parse_ROIs(objects=blue_objects, grp_col='blue', c1=comp1, c2=comp2,
                                                      im=aligned_image.loc[:, frame, :, :],
                                                      annotated_im=annotated_image, stats=stats_table)
        except TypeError:
            print("No blue huh?")


        stats_columns = ["Group", "ROI"]
        stats_columns.extend(aligned_image.species.values)
        stats_columns.append("Ratio_" + args.compare1 + "x" + args.compare2)
        stats_columns = [re.sub(" ","_", x) for x in stats_columns]
        stats_table = pd.DataFrame(stats_table, columns=stats_columns)

        stats_table.to_csv(args.output + "_f0" +
                           "_ratio" + re.sub(" ", "_", args.compare1) +"-x-" + re.sub(" ", "_", args.compare2) +
                           ".tsv", sep="\t", index=False)
        plt.clf()
        plt.imshow(annotated_image, interpolation='none', cmap='gray')
        plt.axis('off')
        plt.colorbar().set_label("Fraction (per ROI)")
        plt.title(sims.utils.format_species(args.compare1) + " / (" + sims.utils.format_species(args.compare1)+" + "+
                  sims.utils.format_species(args.compare2) + ")")
        plt.savefig(fname=args.output + "_f" + str(frame) +
                          "_ratio" + re.sub(" ", "_", args.compare1) +"-x-"+ re.sub(" ", "_", args.compare2) + ".png")
        if not args.no_show:
            plt.show()


