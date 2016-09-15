
Comparison of Multivariate Gaussian transformations
---------------------------------------------------

Example code for presenting the comparison of MVN using: sampling, linearization (like EKF) and unscented transformation (like UKF).

The code simplifies the testing for different conditions and provides a 2x3 plot with input in the rows, and outputs in the columns (flattened over 2D).

Metrics for comparing the resulting gaussians are also provided: KL, BC, DB and Frichet. Note that Frichet depends on chol so it is less reliable with semidefinite variances.

Bonus: we have added a metrics of linearity of functions based on the Hessian. This will allow to evaluate the quality of linearization based on the specific non-linearity. The Hessian of a function f:Rn -> Rm is a tensor: [m,n,n].

Example below of Function: [x2; 2x1(x2+x1)] evaluated at the Gaussian with mean [0.0054,1.01] and covariance [0.1525,-0.0823,",0.0575]. Sampling done with 500 points.

![Example of result for the function [x2; 2x1(x2+x1)] with point [0.0054,1.01] and variance [0.1525,-0.0823;-0.0823,0.0575] with sampling of 500](example.png)

Related:
* Merwe PhD Thesis about UKF
* Lorenzo Peppoloni PhD Thesis where discussing this issue for our motion reconstruction works

References
------------------
Please cite the following paper if you use some these consideration of this package for a Scientific Publication. We employed them for Inertial motion reconstruction.

 Peppoloni L., Filippeschi A., Ruffaldi E. & Avizzano C.A. (2013). A novel 7 degrees of freedom model for upper limb kinematic reconstruction based   on wearable sensors. In Intelligent Systems and Informatics (SISY), 2013 IEEE 11th International Symposium on (pp. 105-110). .  doi:10.1109/SISY.2013.6662551
