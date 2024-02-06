from plantcv import plantcv as pcv
import matplotlib
import cv2
import numpy as np
import argparse 
from  matplotlib import pyplot as plt
import os

### Parse command-line arguments
def options():
    parser = argparse.ArgumentParser(description="Imaging processing with opencv")
    parser.add_argument("-r", "--result", help="result file.", required=False)
    parser.add_argument("-i", "--image", help="Input image file.", required=True)
    parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=True)
    parser.add_argument("-n", "--names", help="path to txt file with names of genotypes to split images into", required =False)
    parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action=None)
    args = parser.parse_args()
    return args

def main():
    # Get options
    args = options()

    # Set variables
    pcv.params.debug = args.debug     # Replace the hard-coded debug with the debug flag

    #import image
    img, path, filename = pcv.readimage(args.image)
    
    mask = pcv.naive_bayes_classifier(rgb_img=img, 
                                  pdf_file="/shares/mgehan_share/acasto/brassica/SV_naive_bayes_pdfs.txt")
    
    plant_img = pcv.apply_mask(mask=(mask['plant']), img=img, mask_color='black')
    filled = pcv.fill(bin_img=(mask['plant']), size=800)
    #pcv.print_image(img=plant_img, filename = os.path.join(args.outdir, filename[:-4] + "_maskimg.png"))
    
    id_objects, obj_hierarchy = pcv.find_objects(img=img, mask=filled)
    
    roi_contour, roi_hierarchy = pcv.roi.rectangle(img=img, x=800, y=750, h=400, w=800)
    
    roi_objects, roi_obj_hierarchy, kept_mask, obj_area = pcv.roi_objects(img=img, roi_contour=roi_contour, 
                                                                          roi_hierarchy=roi_hierarchy,
                                                                          object_contour=id_objects,
                                                                          obj_hierarchy=obj_hierarchy, 
                                                                          roi_type='partial')
    
    filtered_contours, filtered_hierarchy, filtered_mask, filtered_area = pcv.roi_objects(
        img=img, roi_type="partial", roi_contour=roi_contour, roi_hierarchy=roi_hierarchy, object_contour=id_objects, 
        obj_hierarchy=obj_hierarchy)
    
    plant_contour, plant_mask = pcv.object_composition(img=img, contours=filtered_contours, hierarchy=filtered_hierarchy)
    
    nir_path = pcv.get_nir(path=path, filename=filename)
    grey, grey_path1, grey_filename = pcv.readimage(filename=nir_path, mode='grey')
    nir, nir_path1, nir_filename = pcv.readimage(filename=nir_path)
    
    nmask = pcv.resize(img=plant_mask, resize_x=0.285, resize_y=0.285)
    #erode_img = pcv.erode(gray_img=nmask, ksize=3, i=1)
    newmask = pcv.crop_position_mask(img=grey, mask=nmask, x=91, y=2, v_pos="top", h_pos="right")
    
    plant_img = pcv.apply_mask(mask=newmask, img=grey, mask_color='white')
    #pcv.print_image(img=plant_img, filename = os.path.join(args.outdir, filename[:-4] + "_maskimg.png"))
    
    nir_objects, nir_hierarchy = pcv.find_objects(img=grey, mask=newmask)
    nir_combined, nir_combinedmask = pcv.object_composition(img=grey, contours=nir_objects, hierarchy=nir_hierarchy)
    nir_hist = pcv.analyze_nir_intensity(gray_img=nir, mask=nir_combinedmask, bins=256, histplot=False)
    
    pcv.print_results(filename=args.result)
    # Clear the measurements stored globally into the Ouptuts class
    pcv.outputs.clear()
    
    
if __name__ == '__main__':
    main()