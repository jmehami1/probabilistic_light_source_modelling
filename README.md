# Probabilistic Light Source Modelling

This repository has a submodule. Make sure to add the `--recurse-submodules` flag when pulling or cloning.

```bash
#clone repo and submodules
git clone --recurse-submodules **REPO**

#first clone repo and then initialise and update submodule
git clone **REPO**
git submodule init
git submodule update
```



## Non-isotropic light source simulation

Relevant script: **main_light_simulation.m**

Simulates a near-field non-isotropic light source using a measured RID and fits the following regression models for its irradiance:

- Least squares
- GP with zero-mean
- GP with constant-mean
- GP with light-source-mean

The results are saved in `sim_results`



## Real light source triangulation

Relevant script: **main_light_triangulation.m**

For a single real light source, determines its principal direction and location (w.r.t to a frame camera) if it was assumed to be a non-isotropic point source.  Need to capture images with the reflective hemisphere attached to a ChArUco board.

 The data directory should have the following structure:

<pre>light_trig
├── Frame
│   ├── img1.png
│   ├── img2.png
│   └── ...
├── Result (Created after running script)
│   └── pt_light.mat</pre>




Results will be saved in a MAT file called **pt_light.mat** located in a new directory called `Result` 

| Variable Name | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| locLightSrc   | [3x1] Translation vector of the point source location w.r.t to the RGB frame camera's coordinate system in world units (metres) |
| rotLightSrc   | [3x3] Rotation matrix of the point source w.r.t to the RGB frame camera's coordinate system. The Z-axis is inline with the principal direction of the source. |



You can edit the corresponding parameter file located in`parameter_files/light_triangulation.yaml` .  *(You should not need to change the board parameters unless the board has been changed)*

| Flags       | Description                                                  |
| ----------- | ------------------------------------------------------------ |
| DisplayOn   | [true/false] Turns the intermediate visualisations on/off. Turn off to speed up processing |
| FrameCamera | The name of the frame camera. This name must match the filename of the frame camera intrinsic .mat file |
| UseTestData | [true/false] Using the test data present in directory **test_data** |



## Real light source modelling

Relevant script: **main_light_modelling.m**

For a single light source, models the emitted spatial irradiance for the following regression models:

- Least squares
- GP with light-source-mean

Need to capture images of the diffuse reflectance stripe attached to an ArUco board. **Light source triangulation must be done prior to running this script.**

 The data directory should have the following structure:

<pre>light_map
├── Frame
│   ├── img1.png
│   ├── img2.png
│   └── ...
├── Line-scan
│   ├── hs1.png
│   ├── hs2.png
│   └── ...
├── Result (Created after running script)
│   ├── gp_lightsrc_optm.mat
│   ├── imgExtrinsicPose.mat
│   └── ...
├── dark_ref.png
├── pt_light.mat (Result from light source triangulation)
└── camera_system_optimised_parameters.mat (Result from camera system calibration)</pre>

**imgExtrinsicPose.mat** saves the detected poses from each image of the board. Requires the script to be run once on the dataset.

The results can be found in the MAT file **gp_lightsrc_optm.mat**  located in a new directory called `Result` . If script is run once, the trained hyper-parameters will be saved and training will not be done in the next run (unless the down-sample rate or STD radiance noise is changed).

| Variable Name    | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| downSamplingGP   | Down-sampling rate applied to training data for current trained parameters |
| hypOptCell       | Structure of trained hyper-parameters and GP training information for each band of the hyperspectral camera (164). |
| stdRadianceNoise | STD of noise in radiance measurements for current trained parameters |



You can edit the corresponding parameter file located in`parameter_files/light_modelling.yaml` .  *(You should not need to change the board parameters unless the board has been changed)*

| Flags                     | Description                                                  |
| ------------------------- | ------------------------------------------------------------ |
| DisplayOn                 | [true/false] Turns the intermediate visualisations on/off. Turn off to speed up processing |
| FrameCamera               | The name of the frame camera. This name must match the filename of the frame camera intrinsic .mat file |
| UseTestData               | [true/false] Using the test data present in directory **test_data** |
| sigmaRadianceNoisePercent | Zero-mean Gaussian STD noise in radiance measurements as a percentage of the measurement. |
| DownSamplingGP            | Down-sampling rate applied to training data for GP regression |





A hyperspectral camera requires a powerful directed light source that will illuminate a surface of interest that is being captured. Understanding how the radiant intensity of this source is distributed across this surface can aid in trying to estimate material properties from the hyperspectral images. This work aims to build a probabilistic point model of a single **near-field non-isotropic** light source using a **Gaussian Process** with a **non-zero mean function**. Non-isotropic behaviour means that the light source has a **principal direction** where radiant intensity is the highest.

The hyperspectral camera we have worked with captures in a line-scan manner, so therefore, we have incorporated an RGB frame camera to compensate for the missing dimension. The frame camera coordinate frame is **always** treated to be the world coordinate frame. This work assumes the camera system has been calibrated prior using calibration scripts from the **Line-scan_Frame_Camera_Calibration** repository.

There are three key components to this work:

1. Realistic light simulator with a real RID curve to validate GP light source model for a single band

2. Real data captured using camera system

   ​	a. Captured frame images to triangulate the location and direction of light source in world coordinates

   ​	b. Captured frame and hyperspectral images to build the GP model w.r.t to the light source.

3. Calculate estimated reflectance from radiance hypercubes

   ​	a. Estimated reflectance using existing Robles and Krebs method with information from GP light source model.

   ​	b. Estimate reflectance using our method of KrebsShape which incorporates shape information (surface normals)



## Light Source Simulator

MATLAB Simulator 

<img align="center" src="Diagrams/light_simulator.png" alt="Light simulator in MATLAB" width="1000" />



## main_light_modelling

This script is used to build the light source model using collected data from the camera system for each of the bands of the hyperspectral camera. Each band will have a GP and least-squares model.

The script can be run using data located in the `test_data` directory. The parameter `UseTestData` will need to be set to **true** in the YAML parameter file located at `\parameter_files\light_source_modelling.yaml`
