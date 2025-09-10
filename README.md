# Simulating radiative transfer using Monte Carlo methods

**Section I**: Random numbers drawn from a non-flat probability distribution, this is achieved once by using a rejection method and again by using an inverse cumulative distribution function. Two histograms are produced to verify that the methods are working accordingly.

**Section II**: Isotropic scattering. 10<sup>6</sup> photons are ejected from the origin at varying angles ϑ, which scatter until they escape the medium. Scattering patterns are plotted. 

**Section III**: Rayleigh scattering. Optical depth is varied to show why the sky is blue.

## Plots

Generated using Maple Software.

### Rejection method (Sec. I)

![Figure 1](/Plots/RejectionMethod.png)

**Figure 1**: Black lines represent y = exp −x. *LEFT*: A scatter plot highlighting data points sampled using the rejection method. The blue data points represent the sampled data. *RIGHT*: A histogram highlighting where a higher density of data is sampled. The bin width is narrow for large y, since more data is generated in this space.
<br/>


### Inverse CDF method (Sec. I)

![Figure 2](/Plots/InverseCDF.png)

**Figure 2**: A histogram displaying data sampled using an inverse CDF. The black line represents y = exp −x. The bins can be seen to follow the trajectory of the black line.
<br/>


### Isotropic scattering (Sec. II)

![Figure 3](/Plots/IsotropicScatterFinalPositions.png)

**Figure 3**: Photons were modelled to scatter isotropically. *LEFT*: Scatter plot of the photons’ final exit positions in the (y, z) plane. *RIGHT*: Scatter plot of the photons’ final exit positions in the (x, y) plane.
<br/>

![Figure 4](/Plots/IsotropicScatterMovement.png)

**Figure 4**: *LEFT*: A scatter plot showing which angles photons likely exit the medium at. Photons were binned according to scattering angle ϑ, the bin width was constant at 0.1 radians. *RIGHT*: A 3D plot showing the scattering patterns of 5 photons from their initial movement from the origin till they escape the medium.


### Rayleigh scattering (Sec. III)
