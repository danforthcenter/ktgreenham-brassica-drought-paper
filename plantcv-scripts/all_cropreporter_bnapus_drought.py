from plantcv import plantcv as pcv
import matplotlib
import cv2
import numpy as np
import argparse
from  matplotlib import pyplot as plt
import os
from skimage.util import img_as_ubyte
from plantcv.parallel import workflow_inputs

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


if (ps.lightadapted is not None) and (ps.darkadapted is not None):
    

    # select the chlorophyll frame
    img = ps.chlorophyll.sel(frame_label = "Chl").data
    pcv.plot_image(img)

    # standard thresholding
    plant_mask = pcv.threshold.binary(gray_img=img, threshold=450, max_value=255, object_type="light")

    # fill in all the little white specks (salt)
    filled_mask = pcv.fill(bin_img=plant_mask, size=500)

    #find objects
    obj, obj_hierarchy = pcv.find_objects(img=img, mask=filled_mask)

    # Specify a list of coordinates of desired ROIs 
    roi_contours = pcv.roi.multi(img=img, coord=[(400,510), (890,510)],radius=150)

    #label rois
    all_masks = []
    labels = ["control", "drought"]

    #loop through measurements on both plants separately
    i=0
    for roi, hierarchy in roi_contours:
        # Find objects
        filtered_contours, filtered_hierarchy, filtered_mask, filtered_area = pcv.roi_objects(
            img=img, roi_type="partial", roi_contour=roi, roi_hierarchy=hierarchy, object_contour=obj, 
            obj_hierarchy=obj_hierarchy)

        # Combine objects together in each plant     
        plant_contour, plant_mask = pcv.object_composition(img=img, contours=filtered_contours, hierarchy=filtered_hierarchy)        

        # Analyze the shape of each plant 
        analysis_images = pcv.analyze_object(img=img, obj=plant_contour, mask=plant_mask, label=labels[i])

        # Save the image with shape characteristics 
        img_copy = analysis_images

        # get dark adapted data
        dark_da, dark_fig, dark_df = pcv.photosynthesis.reassign_frame_labels(ps_da=ps.darkadapted, mask=plant_mask)

        # get light adapted data
        light_da, light_fig, light_df = pcv.photosynthesis.reassign_frame_labels(ps_da=ps.lightadapted, mask=plant_mask)

        # calculate Fv/Fm
        fvfm, fvfm_hist = pcv.photosynthesis.analyze_yii(ps_da=dark_da, mask=plant_mask,measurement_labels=["Fv/Fm"], label=labels[i])

        # makes an Fv/Fm pseudocolor plant image
        fvfm_cmap = pcv.visualize.pseudocolor(gray_img=fvfm.data.squeeze(), mask=plant_mask, cmap="viridis", min_value=0.5, max_value=1, title="Fv/Fm", background = "black")
    
        fvfm_cmap.savefig(f"{imgname}_fvfm_pseudo_{labels[i]}")


        #chlorophyll index
        ci = pcv.spectral_index.ci_rededge(hsi=ps.spectral)
        ci_hist = pcv.hyperspectral.analyze_index(index_array=ci, mask=plant_mask, min_bin=0, max_bin=20, label=labels[i])

        # fq/fm 
        fqfm, fqfm_hist = pcv.photosynthesis.analyze_yii(ps_da=light_da, mask=plant_mask, measurement_labels=["Fq'/Fm'"], label=labels[i])

        fqfm_cmap = pcv.visualize.pseudocolor(gray_img=fqfm.data.squeeze(), mask=plant_mask, cmap="viridis", min_value=0, max_value=1, title="Fq'/Fm'")

        fqfm_cmap.savefig(f"{imgname}_fqfm_pseudo_{labels[i]}")

        # NPQ
        npq, npq_hist = pcv.photosynthesis.analyze_npq(ps_da_light=light_da, ps_da_dark=dark_da, mask=plant_mask, min_bin = 0, max_bin = 5, measurement_labels=["NPQ"], label=labels[i])

        npq_cmap = pcv.visualize.pseudocolor(gray_img=npq.data.squeeze(), mask=plant_mask, cmap="viridis", min_value=0, max_value=1, title="NPQ")

        npq_cmap.savefig(f"{imgname}_npq_pseudo_{labels[i]}")

        # Anthocyanin Reflective Index
        ari = pcv.spectral_index.ari(hsi=ps.spectral)
        ari_ps = pcv.visualize.pseudocolor(gray_img=ari.array_data, min_value=0, max_value=100, cmap="Purples", mask=plant_mask, background="black", title="Anthocyanin Reflectance Index")
        
        ari_ps.savefig(f"{imgname}_ari_pseudo_{labels[i]}")


        ari_hist = pcv.hyperspectral.analyze_index(index_array=ari, mask=plant_mask, min_bin=-10, max_bin=100, label=labels[i])

        # Chlorophyll Index Red Edge
        ci = pcv.spectral_index.ci_rededge(hsi=ps.spectral)
        ci_ps = pcv.visualize.pseudocolor(gray_img=ci.array_data, min_value=0, max_value=10, cmap="Greens", mask=plant_mask, background="black", title="Chlorophyll Index Red Edge")

        ci_ps.savefig(f"{imgname}_ci_pseudo_{labels[i]}")

        ci_hist = pcv.hyperspectral.analyze_index(index_array=ci, mask=plant_mask, min_bin=0, max_bin=15, label=labels[i])

        # Normalized Difference Vegetation Index
        ndvi = pcv.spectral_index.ndvi(hsi=ps.spectral)
        ndvi_ps = pcv.visualize.pseudocolor(gray_img=ndvi.array_data, min_value=0, max_value=1, cmap="jet", mask=plant_mask, background="black", title="Normalized Difference Vegetation Index")

        ndvi_ps.savefig(f"{imgname}_ndvi_pseudo_{labels[i]}")

        ndvi_hist = pcv.hyperspectral.analyze_index(index_array=ndvi, mask=plant_mask, min_bin=0, max_bin=1, label=labels[i])

        i = i+1



# Write shape and PSII data to results file
pcv.outputs.save_results(filename=args.result, outformat="json")

# Clear the measurements stored globally into the Ouptuts class
pcv.outputs.clear()
