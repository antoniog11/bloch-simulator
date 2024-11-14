# 3D Bloch Simulator using Matlab ODE Solver
This code allows you to simulate and visualize the precession behavior of the net magnetization vector that is the cornerstone of NMR and MRI. The primary file is precession.m which contains the user controlled simulation parameters, instantiation of the ODE solver function, and 3D plotting of the magnetization vector and its vector components. 

# Examples
![](/Examples/free-precession-no-relaxation-lab.mp4)

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

