////////////////////////////////////////////////////////////////
// Time-Lapse Fusion implementation
//
// Copyright (C) 2012, Francisco J. Estrada
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Main module implementing Time-Lapse Fusion over a sequence of
// images as described in the CPCV-2012 paper. Based on the Exposure
// Fusion algorithm by Mertens et al. (see ExposureFuse() below).
//
////////////////////////////////////////////////////////////////

#include "imageProc.h"

// Global program data - Basically all the work data needed by the
// fusion function.
char srcPath[2048];			// Source image path
char dstPath[2048];			// Output image path
double alphaC,alphaS,alphaE;		// Alpha values for contrast, saturation
					// and exposedness.
double tau;				// Tau for temporal blending
int maxFrames;				// Max number of frames to keep=round(3*tau)
struct pyramid **srcPyr;		// Array of Laplacian pyramids for source images
struct pyramid **wgtPyr;		// Array of Gaussian pyramids for weight maps
double *weights;			// Array of weights for temporal blending
double sigma;

void ExposureFuse(void)
{
 // Standard Exposure Fusion as described in
 //
 // Exposure Fusion
 // Tom Mertens, Jan Kautz and Frank Van Reeth
 // Pacific Graphics 2007
 //
 // This code allows for incremental computation so arbitrarily long
 // image sequences can be processed without concern for memory
 // limitations (and that's the only reason we have a separate function
 // from TimeLapseFuse(), otherwise we could use TimeLapseFuse()
 // with all temporal decay weights set to 1).
 char s[2048],oname[2048];
 FILE *f;
 struct image *im, *wghts, *t1;
 struct pyramid *pyrI, *pyrW, *tPyr, *fPyr, *wPyr;
 int i;

 // Create a list of source images, whether .ppm or .PPM, and store it
 // in a text file (C is a scripting language!)
 sprintf(&s[0],"find %s -name \"*.ppm\" | sort > inputList.txt",srcPath);
 system(s);
 sprintf(&s[0],"find %s -name \"*.PPM\" | sort >> inputList.txt",srcPath);
 system(s);

 // Open input list for reading and process image sequence
 f=fopen("inputList.txt","r");
 fPyr=NULL;
 wPyr=NULL;

 // Process all input images
 while (fgets(&s[0],2040,f))
 {
  // Must remove the end-of-lines
  for (i=strlen(s)-1;i>0;i--)
   if (s[i]=='m'){s[i+1]='\0';break;}
  fprintf(stderr,"Processing: %s\n",s);

  // Read in new image, process image, weights, and pyramids
  im=readPPM(s);
  if (im==NULL)
  {
   fprintf(stderr,"ExposureFuse(): Error, can not read input image %s\n",s);
   return;
  }
  // Compute weights for this image
  wghts=computeWeightMap(im,alphaC,alphaS,alphaE);
  if (!wghts)
  {
   fprintf(stderr,"ExposureFusion: Error, unable to compute weight map\n");
   deleteImage(im);
   return;
  }

  // pyrI, pirW, and tPyr are temporary pyramids for holding the image's
  // Laplacian pyramid, Gaussian weights pyramid, and weighted Laplacian
  // pyramid respectively.
  // wPyr and fPyr are the accumulator pyramids for the final blended
  // frame. They must be allocated when we process the first image.
  pyrI=LaplacianPyr(im,5);	// Laplacian pyramid for new input image
  if (!pyrI)
  {
   fprintf(stderr,"ExposureFuse(): Error, unable to compute Laplacian pyramid\n");
   deleteImage(im);
   deleteImage(wghts);
   return;
  }
  pyrW=GaussianPyr(wghts,5);	// Gaussian pyramid of weight map
  if (!pyrW)
  {
   fprintf(stderr,"ExposureFuse(): Error, unable to compute Gaussian pyramid\n");
   deleteImage(im);
   deleteImage(wghts);
   deletePyramid(pyrI);
   return;
  }
  tPyr=weightedPyr(pyrI,pyrW);	// Weighted Laplacian pyramid
  if (!tPyr)
  {
   fprintf(stderr,"ExposureFuse(): Error, unable to compute weighted pyramid\n");
   deleteImage(im);
   deleteImage(wghts);
   deletePyramid(pyrI);
   deletePyramid(pyrW);
   return;
  }
  // Un-weighted Laplacian pyramid, image, and weight map are not needed anymore
  deletePyramid(pyrI);
  deleteImage(im);
  deleteImage(wghts);

  if (fPyr==NULL)	// First frame, allocate data for blended pyramid
  {
   fPyr=(struct pyramid *)calloc(1,sizeof(struct pyramid));
   wPyr=(struct pyramid *)calloc(1,sizeof(struct pyramid));
   if (!fPyr||!wPyr)
   {
    fprintf(stderr,"ExposureFuse(): Can not allocate memory for blended pyramid\n");
    deletePyramid(tPyr);
    deletePyramid(pyrW);
    return;
   }
   fPyr->levels=tPyr->levels;
   fPyr->images=(struct image **)calloc(tPyr->levels,sizeof(struct image));
   wPyr->levels=tPyr->levels;
   wPyr->images=(struct image **)calloc(tPyr->levels,sizeof(struct image));
   if (!fPyr->images||!wPyr->images)
   {
    fprintf(stderr,"ExposureFuse(): Can not allocate memory for blended pyramid\n");
    free(fPyr);
    free(wPyr);
    deletePyramid(tPyr);
    deletePyramid(pyrW);
    return;
   }
   for (i=0;i<wPyr->levels;i++)
   {
    *(fPyr->images+i)=newImage((*(tPyr->images+i))->sx,(*(tPyr->images+i))->sy,3);
    *(wPyr->images+i)=newImage((*(tPyr->images+i))->sx,(*(tPyr->images+i))->sy,1);
    if (!(*(wPyr->images+i))||!(*(fPyr->images+i)))
    {
     fprintf(stderr,"ExposureFuse(): Unable to allocate image data for blended pyramid\n");
     free(fPyr->images);
     free(wPyr->images);
     free(fPyr);
     free(wPyr);
     deletePyramid(tPyr);
     deletePyramid(pyrW);
     return;
    }
   }
  }		// First frame data allocated!

  // Accumulate data onto blended frame pyramid and weights pyramid
  for (i=0;i<tPyr->levels;i++)
  {
   pointwise_add(*(fPyr->images+i),*(tPyr->images+i));		// Weighted Laplacian
   pointwise_add(*(wPyr->images+i),*(pyrW->images+i));		// Weight map
  }

  // Clean up for this iteration.
  deletePyramid(tPyr);
  deletePyramid(pyrW);
 } // End while
 fclose(f);		// Done reading from input image list

 // Divide by weights-sum map to obtain final frame
 for (i=0;i<fPyr->levels;i++)
 {
  t1=newImage((*(fPyr->images+i))->sx,(*(fPyr->images+i))->sy,3);
  memcpy(t1->layers[0],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
  memcpy(t1->layers[1],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
  memcpy(t1->layers[2],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
  pointwise_div(*(fPyr->images+i),t1);
  deleteImage(t1);
 }
 // Reconstruct!
 t1=collapsePyr(fPyr);

 // Output blended frame
 sprintf(&oname[0],"%sExpFusionOutput.ppm",dstPath);
 writePPM(oname,t1);
 deleteImage(t1);

 // Final clean-up
 deletePyramid(fPyr);
 deletePyramid(wPyr);
}

void TimeLapseFuse(void)
{
 // Perform time-lapse fusion on the input sequence.
 char s[2048],oname[2048];
 FILE *f;
 struct image *im, *wghts, *t1;
 struct pyramid *pyrI, *pyrW, *tPyr, *fPyr, *wPyr;
 int i,j,frameNo;

 // Create a list of source images, whether .ppm or .PPM, and store it
 // in a text file (C is a scripting language!)
 sprintf(&s[0],"find %s -name \"*.ppm\" | sort > inputList.txt",srcPath);
 system(s);
 sprintf(&s[0],"find %s -name \"*.PPM\" | sort >> inputList.txt",srcPath);
 system(s);

 // Compute weights for image blending depending on frame age
 // Weight pyramids are progressively scaled to achieve the
 // desired blending weight.
 for (i=0;i<maxFrames;i++)
  *(weights+i)=exp(-(i*i)/(2*sigma*sigma));

 // Open input list for reading and process image sequence
 frameNo=0;
 f=fopen("inputList.txt","r");
 while (fgets(&s[0],2040,f))
 {
  // Must remove the end-of-lines
  for (i=strlen(s)-1;i>0;i--)
   if (s[i]=='m'){s[i+1]='\0';break;}
  fprintf(stderr,"Processing: %s\n",s);

  // Delete last (oldest) image and associated data if it exists
  if (*(srcPyr+maxFrames-1)) deletePyramid(*(srcPyr+maxFrames-1));
  if (*(wgtPyr+maxFrames-1)) deletePyramid(*(wgtPyr+maxFrames-1));

  // Shift remaining images/weights/pyramids down by 1
  for (i=maxFrames-1;i>0;i--)
  {
   *(srcPyr+i)=*(srcPyr+i-1);
   *(wgtPyr+i)=*(wgtPyr+i-1);
  }

  // Read in new image, process image, weights, and pyramids
  im=readPPM(s);
  if (im==NULL)
  {
   fprintf(stderr,"TimeLapseFuse(): Error, can not read input image %s\n",s);
   return;
  }

  wghts=computeWeightMap(im,alphaC,alphaS,alphaE);
  if (!wghts)
  {
   fprintf(stderr,"TimeLapseFusion: Error, unable to compute weight map\n");
   return;
  }

  pyrI=LaplacianPyr(im,5);	// Laplacian pyramid for new input image
  if (!pyrI)
  {
   fprintf(stderr,"TimeLapseFusion: Error, unable to compute Laplacian pyramid\n");
   return;
  }
  deleteImage(im);

  pyrW=GaussianPyr(wghts,5);	// Gaussian pyramid of weight map
  if (!pyrW)
  {
   fprintf(stderr,"TimeLapseFusion: Error, unable to compute Gaussian pyramid\n");
   deletePyramid(*(srcPyr));
   return;
  }
  *(wgtPyr)=pyrW;	// Gaussian pyramid at first entry in wgtPyr
  deleteImage(wghts);

  tPyr=weightedPyr(pyrI,pyrW);	// Weighted Laplacian pyramid
  if (!tPyr)
  {
   fprintf(stderr,"TimeLapseFusion: Error, unable to compute weighted pyramid\n");
   deletePyramid(pyrW);
   deletePyramid(pyrI);
   return;
  }
  *(srcPyr)=tPyr;	// Weighted Laplacian pyramid for image at first entry in srcPyr
  deletePyramid(pyrI);	// Un-weighted Laplacian pyramid is not needed!

  // Create empty pyramids for blending and blending weights
  fPyr=(struct pyramid *)calloc(1,sizeof(struct pyramid));
  wPyr=(struct pyramid *)calloc(1,sizeof(struct pyramid));
  if (!fPyr||!wPyr)
  {
   fprintf(stderr,"TimeLapseFuse(): Out of memory!\n");
   fprintf(stderr,"Bailing out... this will leak memory!\n");
   return;
  }
  fPyr->levels=tPyr->levels;
  wPyr->levels=tPyr->levels;
  fPyr->images=(struct image **)calloc(tPyr->levels,sizeof(struct image *));
  wPyr->images=(struct image **)calloc(tPyr->levels,sizeof(struct image *));
  if (!fPyr->images||!wPyr->images)
  {
   fprintf(stderr,"TimeLapseFuse(): Out of memory!\n");
   fprintf(stderr,"Bailing out... this will leak memory!\n");
   return;
  }
  for (i=0;i<tPyr->levels;i++)
  {
   *(fPyr->images+i)=newImage((*(tPyr->images+i))->sx,(*(tPyr->images+i))->sy,3);
   *(wPyr->images+i)=newImage((*(tPyr->images+i))->sx,(*(tPyr->images+i))->sy,1);
   if (!(*(fPyr->images+i))||!(*(wPyr->images+i)))
   {
    fprintf(stderr,"TimeLapseFuse(): Out of memory!\n");
    fprintf(stderr,"Bailing out... this will leak memory!\n");
    return;
   }
  }

  // Blend!
  for (i=0;i<maxFrames;i++)			// For each of the images in current window
  {
   // Accumulate weighted Laplacian and weight map onto blended result pyramids
   tPyr=*(srcPyr+i);
   if (!tPyr) break;
   for (j=0;j<fPyr->levels;j++)
   {
    t1=copyImage(*(tPyr->images+j));
    image_scale(t1,*(weights+i));
    pointwise_add(*(fPyr->images+j),t1);
    deleteImage(t1);
   }
   tPyr=*(wgtPyr+i);
   for (j=0;j<wPyr->levels;j++)
   {
    t1=copyImage(*(tPyr->images+j));
    image_scale(t1,*(weights+i));
    pointwise_add(*(wPyr->images+j),t1);
    deleteImage(t1);
   }
  }	// End for
  for (i=0; i<fPyr->levels; i++)	// Divide by accumulated weight map at each level
  {
   t1=newImage((*(wPyr->images+i))->sx,(*(wPyr->images+i))->sy,3);
   memcpy(t1->layers[0],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
   memcpy(t1->layers[1],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
   memcpy(t1->layers[2],(*(wPyr->images+i))->layers[0],t1->sx*t1->sy*sizeof(double));
   pointwise_div(*(fPyr->images+i),t1);
   deleteImage(t1);
  }

  // Done, reconstruct final blended frame
  t1=collapsePyr(fPyr);
  // Output blended frame
  sprintf(&oname[0],"%sTLF_%09d.ppm",dstPath,frameNo);
  writePPM(oname,t1);
  frameNo++;
  // Clean-up
  deleteImage(t1);
  deletePyramid(fPyr);
  deletePyramid(wPyr);
 } // End while

 fclose(f);		// Done reading from input image list

 // Clean-up after ourselves!
 for (i=0;i<maxFrames;i++)
 {
  if (*(srcPyr+i)) deletePyramid(*(srcPyr+i));
  if (*(wgtPyr+i)) deletePyramid(*(wgtPyr+i));
 }
 // All done!
 return;
}

// Main is really only parsing input parameters and doing some cleanup.
// All the actual work is done in either of two functions:
// ExposureFuse()
// or
// TimeLapseFuse()
// Which get called when tau=-1, or tau>=0 respectively.
int main(int argc, char *argv[])
{
 // Main Time-Lapse Fusion module - input parsing and some memory cleanup at the end.

 if (argc!=7)
 {
  fprintf(stderr,"TimeLapseFusion - Wrong number of parameters.\n");
  fprintf(stderr,"Usage: TimeLapseFusion src_path dest_path alphaC alphaS alphaE tau\n");
  fprintf(stderr,"          src_path: Path to source (.ppm) images\n");
  fprintf(stderr,"          dst_path: Output image path\n");
  fprintf(stderr,"          alphaC: Alpha for contribution of the contract map\n");
  fprintf(stderr,"          alphaS: Alpha for contribution of the saturation map\n");
  fprintf(stderr,"          alphaE: Alpha for contribution of the well-exposedness map\n");
  fprintf(stderr,"          tau: Number of frames to blend into each output image\n");
  fprintf(stderr,"Example: TimeLapseFusion ./src/  ./output/ 1.0 1.0 .15 15\n");
  fprintf(stderr,"          if tau=-1, no temporal blending is done, which\n");
  fprintf(stderr,"          results in standard exposure fusion on the input frames\n");
  exit(0);
 }

 strncpy(&srcPath[0],argv[1],2040);
 strncpy(&dstPath[0],argv[2],2040);
 alphaC=atof(argv[3]);
 alphaS=atof(argv[4]);
 alphaE=atof(argv[5]);
 tau=atof(argv[6]);

 // Perform some simple validation on input parameters
 if (alphaC<0||alphaS<0||alphaE<0||alphaC>10||alphaS>10||alphaE>10)
 {
  fprintf(stderr,"Error: Alpha values must be in [0,10]\n");
  exit(0);
 }

 if (tau>=0)
 {
  maxFrames=(int)round(tau);
  sigma=tau/3;
 }
 else
 {
  maxFrames=-1;
  sigma=-1;
 }

 fprintf(stderr,"Time-Lapse Fusion called with\n");
 fprintf(stderr,"Source image directory: %s\n",srcPath);
 fprintf(stderr,"Output image directory: %s\n",dstPath);
 fprintf(stderr,"Alphas [C,S,E]=[%f,%f,%f]\n",alphaC,alphaS,alphaE);
 fprintf(stderr,"Input frames blended into each output frame (-1 is all)=%d\n",maxFrames);
 fprintf(stderr,"Sigma=%f\n",sigma);

 // Initialize arrays
 if (maxFrames>0)
 {
  srcPyr=(struct pyramid **)calloc(maxFrames,sizeof(struct pyramid *));
  wgtPyr=(struct pyramid **)calloc(maxFrames,sizeof(struct pyramid *));
  weights=(double *)calloc(maxFrames,sizeof(double));
  if (!srcPyr||!wgtPyr||!weights)
  {
   fprintf(stderr,"TimeLapseFusion - main(): Out of memory!\n");
   exit(0);
  }
 }

 if (tau>=0) TimeLapseFuse();
 else ExposureFuse();

 // Release memory allocated to arrays - Note that array contents must
 // be de-allocated prior to this!
 if (maxFrames>0)
 {
  free(srcPyr);
  free(wgtPyr);
  free(weights);
 }

 // Done!
 exit(0);
}
