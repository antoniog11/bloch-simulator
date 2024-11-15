# 3D Bloch Simulator Educational Tool for MRI
This code is intended to be used as an educational tool for understanding magnetic resonance imaging. It allows you to play with and visualize the precession behavior of the net magnetization vector of atomic nuclei which is the origin of the signal in NMR and MRI. MRI can be entirely understood from the classical behavior of a magnetic vector, $\vec{M}$ rotating in 3D space under the influence of an external magnetic field $\vec{B_0}$, much like that of a gyroscope rotating in the presence of a graviational field. This simulator can help you visualize and get an intuition for precession, relaxation, and the effect that radio frequency pulses, $\vec{B_1}$, have on the magnetic vector.

The primary file is ```precession.m``` which contains the user controlled simulation parameters at the top, instantiation of the Matlab ODE solver function, and 3D plotting of the magnetization vector (see Examples Below) and aloing with its vector components. 

# Examples

### Free Precession, no relaxation, viewed in laboratory (non-rotating) reference frame.
https://github.com/user-attachments/assets/76b4726d-9115-4f31-8a42-5d03655e25ce

### Free Precession, no relaxation, viewed in rotating reference frame.
https://github.com/user-attachments/assets/87728cba-6b5a-4364-b30c-95fc62755c2e

### Free Precession, with relaxation, laboratory frame
https://github.com/user-attachments/assets/7d96023f-b39e-43db-ab28-b3e0531f5181

### Free Precession, with relaxation, rotating frame
https://github.com/user-attachments/assets/dda76ad1-7a4c-4bf3-b8ef-adc39fe48730

### Forced Precession (applied $\vec{B_1}$), no relaxation, laboratory frame
https://github.com/user-attachments/assets/bb9734ce-e661-40d8-81f9-b6994bd435da

### Forced Precession, no relaxation, rotating frame
https://github.com/user-attachments/assets/28111938-e27d-4447-a5f5-19af3ecae4c1

### Magnetization components ($x,y,z$) versus time, four $\frac{\pi}{2}$ pulses
![](Examples/Example-component-plot.png)

User's can control parameters such as:
- T1 & T2.
- Number of $\vec{B_1}$ excitations.
- $\vec{B_1}$ Magnitude.
- Flip Angle.
- Whether to view the bulk magnetization vector in the rotating or laboratory frame.
- Whether to simulate relaxation or not.
- Whether to trace the tip of the magnetization vector as it moves in 3D space.
- And more!

# Dependencies
- [arrow3d](https://www.mathworks.com/matlabcentral/fileexchange/71994-arrow-3d) on Matlab File exchange to draw the spin vector in 3D space.

# References
This code also makes use of some some code (```showspins.m```) from the [RAD229](https://web.stanford.edu/class/rad229/index.html) course at Stanford taught by Brian Hargreaves and Daniel Ennis.

