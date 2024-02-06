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
    img, path, img_filename = pcv.readimage(args.image)
    filename = img_filename
    
    mask = pcv.naive_bayes_classifier(rgb_img=img, 
                                  pdf_file="/shares/mgehan_share/acasto/brassica/SV_naive_bayes_pdfs.txt")
    
    #plant_img = pcv.apply_mask(mask=(mask['plant']), img=img, mask_color='black')
    #pcv.print_image(img=plant_img, filename = os.path.join(args.outdir, filename[:-4] + "_maskimg.png"))
    
    id_objects, obj_hierarchy = pcv.find_objects(img=img, mask=mask['plant'])
    
    roi_contour, roi_hierarchy = pcv.roi.rectangle(img=img, x=500, y=50, h=1050, w=1000)
    
    roi_objects, roi_obj_hierarchy, kept_mask, obj_area = pcv.roi_objects(img=img, roi_contour=roi_contour, 
                                                                          roi_hierarchy=roi_hierarchy,
                                                                          object_contour=id_objects,
                                                                          obj_hierarchy=obj_hierarchy, 
                                                                          roi_type='partial')
    
    filtered_contours, filtered_hierarchy, filtered_mask, filtered_area = pcv.roi_objects(
        img=img, roi_type="partial", roi_contour=roi_contour, roi_hierarchy=roi_hierarchy, object_contour=id_objects, 
        obj_hierarchy=obj_hierarchy)
    
    plant_contour, plant_mask = pcv.object_composition(img=img, contours=filtered_contours, hierarchy=filtered_hierarchy)
    
    analysis_images = pcv.analyze_object(img=img, obj=plant_contour, mask=plant_mask)
    analysis_image_color = pcv.analyze_color(img, plant_mask, 'hsv')
    
    pcv.print_results(filename=args.result)
    # Clear the measurements stored globally into the Ouptuts class
    pcv.outputs.clear()
    
    
if __name__ == '__main__':
    main()