#include <string>
#include <iostream>
#include <opencv2/aruco.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "mex.h"
#include <vector>
#include <map>

#include "helper_functions.h"

#define NUM_INPUTS_MAT 1
#define NUM_OUTPUTS_MAT 2


/*
 * Used to detect the pixel locations for each corner of each ArUco marker.
*/
cv::Mat DetectArucoMarkerPixel(cv::Mat &image, std::vector<int> &markerIds, std::vector<std::vector<cv::Point2f> > &markerCorners)
{
    cv::Mat imageCopy;
    image.copyTo(imageCopy);

    cv::Ptr<cv::aruco::Dictionary> dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_4X4_50);
    cv::Ptr<cv::aruco::DetectorParameters> params = cv::aruco::DetectorParameters::create();
    params->cornerRefinementMethod = cv::aruco::CORNER_REFINE_SUBPIX;

    cv::aruco::detectMarkers(image, dictionary, markerCorners, markerIds, params);

    if (markerIds.size() > 0)
        cv::aruco::drawDetectedMarkers(imageCopy, markerCorners, markerIds);

    return imageCopy;
}


/*
 * Mex function which interfaces the ArUco board code from opencv to allow it to be run in MATLAB.
 * Inputs prhs[]:
 * [0] - uint8 input image (rows x cols x channels)
 *
 * Outputs plhs[]:
 * [0] - mxDouble identified marker IDs [N x 1]
 * [1] - mxDouble identified marker corner pixel locations [N x 8]
 * [2] - uint8 output image with identified markers (rows x cols x channels) [OPTIONAL]
*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // Check for proper number of input arguments
    if (nrhs != NUM_INPUTS_MAT )
        mexErrMsgTxt("incorrect input arguments");


    //check for proper number of output arguments
    if (nlhs != NUM_OUTPUTS_MAT && nlhs != NUM_OUTPUTS_MAT+1)
        mexErrMsgTxt("incorrect output arguments");

    //**************1ST INPUT**************************

    cv::Mat inImgMat = mx_Array_Image2_Mat(prhs[0]);

    const mwSize *inImgDim = mxGetDimensions(prhs[0]);
    mwSize numImgDim = mxGetNumberOfDimensions(prhs[0]);

    //convert dimensions to integer
    const int inImgH = inImgDim[0];
    const int inImgW = inImgDim[1];

    if (inImgMat.empty())
        mexErrMsgTxt("Could not read in image to opencv MAT type");

    std::vector<int> markerIds;
    std::vector<std::vector<cv::Point2f> > markerCorners;

    //detect all possible markers
    cv::Mat imgMarkerDet = DetectArucoMarkerPixel(inImgMat, markerIds, markerCorners);


    //**************1ST OUTPUT**************************

    int rows = markerIds.size();

    plhs[0] = mxCreateNumericMatrix(rows, 1, mxINT8_CLASS, mxREAL);

    mxInt8* pr = mxGetInt8s(plhs[0]);

    for (int i = 0; i < rows; i++)
        *(pr+i) = markerIds[i];

    //**************SECOND OUTPUT**************************

    rows = markerCorners.size();

    plhs[1] = mxCreateDoubleMatrix(rows, 8, mxREAL);

    double* pr1 = mxGetPr(plhs[1]);

    int arrIndex = 0;

    for (int i = 0; i < 4; i++){
        arrIndex = 0;

        for (int j = 0; j < rows; j++){
            cv::Point2f curCorner = markerCorners[j][i];

            pr1[arrIndex + 2*i*rows] = curCorner.x;
            pr1[arrIndex + (2*i +1)*rows] = curCorner.y;
            arrIndex++;
        }
    }

    //**************THIRD OUTPUT************************** [OPTIONAL]

    if (nlhs == NUM_OUTPUTS_MAT+1) {
        plhs[2] = mxCreateNumericArray(numImgDim,inImgDim, mxUINT8_CLASS, mxREAL);

        char* outMat = (char*) mxGetData(plhs[2]);

        // grayscale image
        if (numImgDim == 2){
            arrIndex = 0;

            //Store image pixel channel colours into a 1D array used for passing to matlab
            for (int j = 0; j < inImgW; j++){
                for (int i = 0; i < inImgH; i++){
                    outMat[arrIndex] = imgMarkerDet.at<char>(i,j);

                    arrIndex++;
                }
            }
        }
        //RGB image
        else {
            cv::Vec3b pixel;
            arrIndex = 0;

            //Store image pixel channel colours into a 1D array used for passing to matlab
            for (int j = 0; j < inImgW; j++){
                for (int i = 0; i < inImgH; i++){
                    pixel = imgMarkerDet.at<cv::Vec3b>(i,j);

                    outMat[arrIndex] = pixel[2];   //R
                    outMat[inImgH*inImgW+arrIndex] = pixel[1]; //G
                    outMat[2*inImgH*inImgW+arrIndex] = pixel[0]; //B

                    arrIndex++;
                }
            }
        }
    }
}
