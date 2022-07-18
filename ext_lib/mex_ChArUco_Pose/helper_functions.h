#include <string>
#include <iostream>
#include <opencv2/aruco/charuco.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "mex.h"
#include <vector>
#include <map>


//Helper function for used in mx_Array_Image2_Mat
mwIndex subs(const mxArray *inputMatrix, const std::vector<mwIndex>& si)
{
    std::vector<mwIndex> v(si);
    return mxCalcSingleSubscript(inputMatrix, si.size(), (!v.empty() ? &v[0] : NULL));
}

cv::Mat mx_Array_Image2_Mat(const mxArray *inputMatrix)
{

//    uint8_T *inImgArr = mxGetUint8s(inputMatrix);
    const mwSize *dims = mxGetDimensions(inputMatrix);
    mwSize ndims = mxGetNumberOfDimensions(inputMatrix);
//    const mxClassID classID = mxGetClassID(inputMatrix);

    // Create cv::Mat object (of the specified depth), equivalent to mxArray.
    // At this point we create either a 2-dim with 1-channel mat, or a 2-dim
    // with multi-channels mat. Multi-dims case is handled above.
    std::vector<int> d(dims, dims+ndims);
    ndims = (d.size()>2) ? d.size()-1 : d.size();
    const mwSize nchannels = (d.size()>2) ? d.back() : 1;
    int depth = CV_8U;
    std::swap(d[0], d[1]);
    // Copy each channel from mxArray to Mat (converting to specified depth),
    // as in: channels[i] <- cast_to_mat_depth(p_(:,:,i))
    std::vector<cv::Mat> channels(nchannels);
    std::vector<mwSize> si(d.size(), 0);                 // subscript index
    const int type = CV_MAKETYPE(depth, 1); // Source type

    for (mwIndex i = 0; i<nchannels; ++i) {
        si[si.size() - 1] = i;                   // last dim is a channel idx
        void *pd = reinterpret_cast<void*>(
                    reinterpret_cast<size_t>(mxGetData(inputMatrix)) +
                    mxGetElementSize(inputMatrix)*subs(inputMatrix ,si));      // ptr to i-th channel data
        const cv::Mat m(ndims, &d[0], type, pd); // only creates Mat headers
        // Read from mxArray through m, writing into channels[i]
        // (Note that saturate_cast<> is applied, so values are clipped
        // rather than wrap-around in a two's complement sense. In
        // floating-point to integer conversion, numbers are first rounded
        // to nearest integer then clamped).
        m.convertTo(channels[i], CV_MAKETYPE(depth, 1));
        // transpose cv::Mat if needed. We do this inside the loop on each 2d
        // 1-cn slice to avoid cv::transpose limitation on number of channels
        //if (transpose)
        cv::transpose(channels[i], channels[i]);  // in-place transpose
    }

    cv::Mat mat(ndims, &d[0], CV_MAKETYPE(depth, nchannels));

    // Merge channels back into one cv::Mat array
    cv::merge(channels, mat);

    //if RGB image, convert from RGB to BGR format
    if (nchannels == 3)
        cv::cvtColor(mat, mat, cv::COLOR_RGB2BGR);

    return mat;
}

/*
 * Converts a mxArray containing a double type matrix to an opencv Mat double matrix
*/
cv::Mat double_mxArray_matrix2cv_Mat_matrix(const mxArray *inputMatrix)
{
    double* inDoubleArr = mxGetDoubles(inputMatrix);
    const mwSize *inDim = mxGetDimensions(inputMatrix);

    //Size of array from Matlab
    const int rows = inDim[0];
    const int cols = inDim[1];

    //Create Mat_ intermediate storage
    cv::Mat_<double> M(rows,cols);

    int indexArr = 0;

    //put all data into the M matrix
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            M[i][j] = *(inDoubleArr + indexArr);
            indexArr++;
        }
    }

    //cast cv::_Mat to cv::Mat
    return ((cv::Mat) M);
}
