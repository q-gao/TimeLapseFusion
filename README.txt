Time-Lapse Fusion implementation as per the CPCV 2012 paper.

This code is licensed under GPL v3.0. If you plan to use it
for anything other than your own personal amusement, please
take a moment to read the included license terms.

The code should compile and run without modification on Linux,
though you will need the pthread libraries if you want multi-thread
support.

To compile. Simply run the compile.sh script.

Usage:

TimeLapseFusion src_path dst_path alphaC alphaS alphaE tau

 src_path specifies where the input images are to be found
    The code will use 'find' to locate and process any
    '.ppm' or '.PPM' images found within the source path.
    This includes sub-directories of the source path.

 dst_path is the destination directory where output images
    will be placed (also in .ppm format)

 alphaC, alphaS, alphaE - are the relarive weights of the
    contrast, saturation, and well-exposedness components
    of the pixelwise weights as specified by the exposure
    fusion algorithm.  

    For simple exposure fusion (single output image from all
    input frames), a reasonable setting is 1.0, 1.0, 1.0

    For time-lapse fusion, reduce the weight of the well-
    exposedness term (images are assumed to be correctly
    exposed). E.g. 1.0, 1.0, .1

  tau - For time-lapse fusion this indicates the number of
        input frames to be blended into each output frame.
        If tau=-1, the program does simple exposure fusion
        and creates a single output image from all input
        frames.

 Note that the code expects a properly formed .ppm header.
 Linux image utilities such as the GIMP and imageMagick 
 will produce such .ppm headers, but last time I checked
 Matlab's imwrite() will not. 

 To easily convert a set of images in a different format
 to ppm, use the mogrify command:

 e.g.     mogrify -format ppm *.jpg

 Similarly the output .ppm files can be easily converted
 to .jpg for movie encoding.

For comments or bug reports, please contact the author

festrada [ a.t. ] utsc [ dot ] utoronto [ dot ] ca

Created by F. Estrada, Oct. 04, 2012

