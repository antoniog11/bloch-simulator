# 3D Bloch Simulator using Matlab ODE Solver
This code allows you to simulate and visualize the precession behavior of the net magnetization vector that is the cornerstone of NMR and MRI. The primary file is precession.m which contains the user controlled simulation parameters, instantiation of the ODE solver function, and 3D plotting of the magnetization vector (see Examples Below) and its vector components. 

# Examples

### Free Precession, no relaxation, laboratory frame.
https://github.com/user-attachments/assets/76b4726d-9115-4f31-8a42-5d03655e25ce

### Free Precession, no relaxation, rotating frame
https://github.com/user-attachments/assets/87728cba-6b5a-4364-b30c-95fc62755c2e

### Free Precession, with relaxation, laboratory frame
https://github.com/user-attachments/assets/7d96023f-b39e-43db-ab28-b3e0531f5181

### Free Precession, with relaxation, rotating frame
https://github.com/user-attachments/assets/dda76ad1-7a4c-4bf3-b8ef-adc39fe48730

### Forced Precession, no relaxation, laboratory frame
https://github.com/user-attachments/assets/bb9734ce-e661-40d8-81f9-b6994bd435da

### Forced Precession, no relaxation, rotating frame
https://github.com/user-attachments/assets/28111938-e27d-4447-a5f5-19af3ecae4c1

### Magnetization components (x,y,z) versus time, four $\frac{\pi}{2}$ pulses
![](Examples/Example-component-plot.png)

User's can control things such as:
- T1 & T2.
- Number of B1 excitations.
- B1 Magnitude.
- Flip Angle.
- Whether to view the bulk magnetization vector in the rotating or laboratory frame.
- Wheter to simulate relaxation or not.
- Whether to trace the tip of the spin vector as it moves in 3D space.

# Dependencies
- [arrow3d](https://www.mathworks.com/matlabcentral/fileexchange/71994-arrow-3d) on Matlab File exchange to draw the spin vector in 3D space.

# References
I made use of some code (showspins.m) from the [RAD229](https://web.stanford.edu/class/rad229/index.html) course at Stanford taught by Brian Hargreaves and Daniel Ennis.

