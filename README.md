# light_source_modelling
Probabilistic modelling of a near-field non-isotropic light source.



MATLAB Simulator 

<img align="center" src="Diagrams/light_simulator.png" alt="Light simulator in MATLAB" width="1000" />



## main_light_modelling

This script is used to build the light source model using collected data from the camera system for each of the bands of the hyperspectral camera. Each band will have a GP and least-squares model.

The script can be run using data located in the `test_data` directory. The parameter `UseTestData` will need to be set to **true** in the YAML parameter file located at `\parameter_files\light_source_modelling.yaml`



# Completed

- Simulator visualisation is working with point source
- Sample the space with a small target
- Obtain the symmetrical sample about the light source position and direction on the view-plane of the line-scan camera
- Estimate light parameters through LS with both points and symmetric points
- Build the light source model using zero-mean GP with SE kernel without symmetric points
- Build the light source model using zero-mean GP with SE kernel with symmetric points
- Build the light source model using constant-mean GP with SE kernel with symmetric points



# To-Do

- Clean up code
- Write mean function for the light source with light source parameters as hyper-parameters which will need to be trained.
- Build the light source model using new mean function GP with SE kernel with/without symmetric points