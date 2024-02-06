from plantcv import plantcv as pcv
import matplotlib
import cv2
import numpy as np
import argparse
from  matplotlib import pyplot as plt
import os
from skimage.util import img_as_ubyte
from plantcv.parallel import workflow_inputs

pcv.params.dpi = 600

### Parse command-line arguments

# Get options
args = workflow_inputs()

# Set variables
pcv.params.debug = args.debug     # Replace the hard-coded debug with the debug flag


#read in the .INF file 
#this reads in all the types of images you took
imgpath = args.psii
ps = pcv.photosynthesis.read_cropreporter(filename=imgpath)
imgname = os.path.splitext(os.path.split(imgpath)[-1])[0]


#*** I highly recommend taking a chlorophyll picture, only because they are useful for masking ***
# the code below selects the chlorophyll frame

img = ps.chlorophyll.sel(frame_label = "Chl").data
pcv.plot_image(img)

# standard thresholding
plant_mask = pcv.threshold.binary(gray_img=img, threshold=450, max_value=255, object_type="light")

# fill in all the little white specks (salt)
filled_mask = pcv.fill(bin_img=plant_mask, size=500)

obj, obj_hierarchy = pcv.find_objects(img=img, mask=filled_mask)

# Specify a list of coordinates of desired ROIs 
#roi_contours = pcv.roi.multi(img=img, coord=[(400,510), (890,510)],radius=150)
roi_contour, roi_hierarchy = pcv.roi.rectangle(img=img, x=250, y=300, h=400, w=800)


roi_objects, hierarchy, kept_mask, obj_area = pcv.roi_objects(img, roi_contour, roi_hierarchy, obj, obj_hierarchy, 'partial')

# Combine objects together in each plant     
plant_contour, plant_mask = pcv.object_composition(img=img, contours=roi_objects, hierarchy=hierarchy)   

# get dark adapted data
dark_da, dark_fig, dark_df = pcv.photosynthesis.reassign_frame_labels(ps_da=ps.darkadapted, mask=plant_mask)

# get light adapted data
light_da, light_fig, light_df = pcv.photosynthesis.reassign_frame_labels(ps_da=ps.lightadapted, mask=plant_mask)

# calculate Fv/Fm
fvfm, fvfm_hist = pcv.photosynthesis.analyze_yii(ps_da=dark_da, mask=plant_mask, measurement_labels=["Fv/Fm"])

# makes an Fv/Fm pseudocolor plant image
fvfm_cmap = pcv.visualize.pseudocolor(gray_img=fvfm.data.squeeze(), mask=plant_mask, cmap="RdYlBu", min_value=0.65, max_value=0.9, title="Fv/Fm", background = "black")

# fq/fm 
fqfm, fqfm_hist = pcv.photosynthesis.analyze_yii(ps_da=light_da, mask=plant_mask, measurement_labels=["Fq'/Fm'"])

fqfm_cmap = pcv.visualize.pseudocolor(gray_img=fqfm.data.squeeze(), mask=plant_mask, cmap="RdYlBu", min_value=0, max_value=0.75, title="Fq'/Fm'", background = "black")

# NPQ
npq, npq_hist = pcv.photosynthesis.analyze_npq(ps_da_light=light_da, ps_da_dark=dark_da, mask=plant_mask, min_bin = 0, max_bin = 5, measurement_labels=["NPQ"])

npq_cmap = pcv.visualize.pseudocolor(gray_img=npq.data.squeeze(), mask=plant_mask, cmap="jet", min_value=0, max_value=2.5, title="NPQ", background = "black")

# Normalized Difference Vegetation Index
ndvi = pcv.spectral_index.ndvi(hsi=ps.spectral)

ndvi_ps = pcv.visualize.pseudocolor(gray_img=ndvi.array_data, min_value=0, max_value=0.8, cmap="viridis", mask=plant_mask, background="black", title="Normalized Difference Vegetation Index")

ndvi_hist = pcv.hyperspectral.analyze_index(index_array=ndvi, mask=plant_mask, min_bin=0, max_bin=1)


# Chlorophyll Index Red Edge
ci = pcv.spectral_index.ci_rededge(hsi=ps.spectral)
ci_ps = pcv.visualize.pseudocolor(gray_img=ci.array_data, min_value=1, max_value=4, cmap="YlGn", mask=plant_mask, background="black", title="Chlorophyll Index Red Edge")

ci_hist = pcv.hyperspectral.analyze_index(index_array=ci, mask=plant_mask, min_bin=0, max_bin=15)

# Anthocyanin Reflective Index
ari = pcv.spectral_index.ari(hsi=ps.spectral)
ari_ps = pcv.visualize.pseudocolor(gray_img=ari.array_data, min_value=-1, max_value=12, cmap="Purples", mask=plant_mask, background="black", title="Anthocyanin Reflectance Index")

ari_hist = pcv.hyperspectral.analyze_index(index_array=ari, mask=plant_mask, min_bin=-10, max_bin=15)

pcv.print_image(fvfm_cmap, f"./bnapus_high_quality_pseudo/{imgname}_fvfm_pseudo.png")

pcv.print_image(fqfm_cmap, f"./bnapus_high_quality_pseudo/{imgname}_fqfm_pseudo.png")

pcv.print_image(npq_cmap, f"./bnapus_high_quality_pseudo/{imgname}_npq_pseudo.png")

pcv.print_image(ndvi_ps, f"./bnapus_high_quality_pseudo/{imgname}_ndvi_pseudo.png")

pcv.print_image(ci_ps, f"./bnapus_high_quality_pseudo/{imgname}_ci_pseudo.png")

pcv.print_image(ari_ps, f"./bnapus_high_quality_pseudo/{imgname}_ari_pseudo.png")


pcv.outputs.save_results(filename=args.result, outformat="json")

# Clear the measurements stored globally into the Ouptuts class
pcv.outputs.clear()