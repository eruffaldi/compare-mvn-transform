
Comparison of Multivariate Gaussian transformations
---------------------------------------------------

Example code for presenting the comparison of MVN using: sampling, linearization (like EKF) and unscented transformation (like UKF).

The code simplifies the testing for different conditions and provides a 2x3 plot with input in the rows, and outputs in the columns (flattened over 2D).

Metrics for comparing the resulting gaussians are also provided: Jeffreys Divergence (aka Symmetric KL), Bhattacharyya Distance, DB, Frichet, Jeffreys Divergence-Riemman, Bhattacharyya-Riemman . Note that Frichet depends on chol so it is less reliable with semidefinite variances. See the paper "A Note on Metric Properties for Some Divergence Measures" ([PDF](http://www.jmlr.org/proceedings/papers/v25/aboumoustafa12/aboumoustafa12.pdf)), Abou-Moustafa and Ferrie,2012 for details.

Bonus: we have added a metrics of linearity of functions based on the Hessian. This will allow to evaluate the quality of linearization based on the specific non-linearity. The Hessian of a function f:Rn -> Rm is a tensor: [m,n,n].

Example below of Function: [x2; 2x1(x2+x1)] evaluated at the Gaussian with mean [0.0054,1.01] and covariance [0.1525,-0.0823,",0.0575]. Sampling done with 500 points.

![Example of result for the function [x2; 2x1(x2+x1)] with point [0.0054,1.01] and variance [0.1525,-0.0823;-0.0823,0.0575] with sampling of 500](CompareNL2.png)

Using Bhattacharyya-Riemman metrics linearization has error 2.73 while unscented 2.5077. Using Symmetric KL divergence linearization error is 0.15549, unscented error is 0.0015486.


![Example of result for the function [exp(x(2)^2); 200*x(1)^2*(x(2)+x(1))] with point [0.0054,1.01] and variance [1E-4,1E-4] with sampling of 500](CompareNL1.png)

Using Bhattacharyya-Riemman metrics linearization has error 5.4201 while unscented 4.8435. Using Symmetric KL divergence linearization is 1.124, unscented is 0.0012436.

Related:
* Merwe PhD Thesis about UKF ([PDF](http://www.cslu.ogi.edu/publications/ps/merwe04.pdf))
* Lorenzo Peppoloni PhD Thesis where discussing this issue for our motion reconstruction works

References
------------------
Please cite the following paper if you use some these consideration of this package for a Scientific Publication. We employed them for Inertial motion reconstruction.

 Peppoloni L., Filippeschi A., Ruffaldi E. & Avizzano C.A. (2013). A novel 7 degrees of freedom model for upper limb kinematic reconstruction based   on wearable sensors. In Intelligent Systems and Informatics (SISY), 2013 IEEE 11th International Symposium on (pp. 105-110). .  doi:10.1109/SISY.2013.6662551
