#### State Estimation/Localization of a stereo camera in 3D using batch Gauss-Newton optimization method
The goal of this project was to estimate the state $(x,y,z,\theta)$ of a vehicle equipped with a stereo camera and inertial measurement unit (IMU) traveling through a field of landmarks. The entire project was completed in MATLAB from a dataset of noisy measurements from the vehicle. 

The batch Gauss-Newton method was implemented to fuse measurements from the IMU with point measurments from the stereo camera. Since the vehicle could freely translate and rotate through the three-dimensional space, the sets of transformation and rotation matricies used to represent the vehicles position over time were matrix Lie groups $SE(3)$ and $SO(3)$. 

This project was completed as a project for one of my graduate courses.
