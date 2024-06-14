# this script can be used to overlay several CT scan slices, saved as 'Figure1.jpg', 'Figure2.jpg.' etc., in the same folder as this script
# the user needs to adapt the number of figures to be merged on line 7 of this script
# after running this script, the merged image CT_scan_merged.jpg will be stored in the same folder as the script and the original images
import cv2

# define the number of figures; this should be at least 2
nr_fig = 3 

figs = {}
for i in range(nr_fig):
    figs[i] = cv2.imread(f'Figure{i+1}.jpg')

merged_image = cv2.addWeighted(figs[0], 0.5, figs[1], 0.5, 0)
if nr_fig > 2:
    for i in range(2, nr_fig):
        merged_image = cv2.addWeighted(merged_image, i / (i + 1), figs[i], 1 / (i + 1), 0)

cv2.imwrite('CT_scan_merged.jpg', merged_image)