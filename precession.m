addpath '/Users/antonioglenn/Desktop/code-and-projects/RAD229 MRI Signals and Sequences/Matlab'
% Bloch Simulator

%% Free Precession in Laboratory Frame without Relaxation (Matrix Operators)
% Animation parameters.
T = 10; % Animation duration (s).
frameRate = 60; % # frames/second.
nFrames = T * frameRate; % Total frames

% Spin parameters
M_0 = [0; 0; 0];              % Origin of origin.
M = [1; 0; 0] + M_0;          % Direciton of Vector w.r.t its origin.
f0 = 1;                       % Larmor frequency (Hz);
dtheta = (2*pi*f0)/frameRate; % Amount of phase per frame [rad]/[frame]). 
scale = 1;

% Plot and appearance.
spinColors = [1 0 0]; 
showspins(M,scale,M_0,spinColors);
xyPlane = [ [-1 1 1 -1]; [-1 -1 1 1]; [-1 -1 -1 -1]];
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);

view(3)
title('3D View')
grid on
axis vis3d % Disable "stretch-to-fill" as view moves

%fill3([-1 1 1 -1], [-1 -1 1 1], [-1 -1 -1 -1], [0.8, 0.5, 0.3],
%'FaceAlpha', 0.3); % Fill xy plane at z = -1
% set(gcf, 'Color', 'black');
% axisColor = 'black';
% axes = plotm(V);
% for i = 1:length(axes)
%     set(axes(i), ...
%         'Color', axisColor, ...
%         'XColor', 'white', ...
%         'YColor', 'white', ...
%         'ZColor', 'white')
% end

% Plot and appearance.
spinColors = [1 0 0]; 
showspins(M,scale,M_0,spinColors);
xyPlane = [ [-1 1 1 -1]; [-1 -1 1 1]; [-1 -1 -1 -1]];
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);

view(3)
title('3D View')
grid on

% Animate precession
for i=1:nFrames
    cla
    M = zrot(rad2deg(dtheta))* M;
    showspins(M,scale,M_0, spinColors);
    pause(0.01);
    fprintf('time = %fs\n', i/frameRate)
end


%% Free Precession in Laboratory Frame without Relaxation (Numerical Solver)

% Animation parameters.
Tend = 5;                   % Animation/Simulation duration (s).
frameRate = 60;             % # frames/second.
time = 0:1/frameRate:Tend;  % Animation/Simulation time points.
nFrames = Tend * frameRate; % Total frames.
gamma = 1;                  % Gyromagnetic ratio for H ([rad]/[s][T]).
B0 = 2*pi;                     % Tesla.
f = (gamma*B0)/(2*pi);      % Larmor frequency (Hz).
M_initial = [0; -0.5; 0.5]; % Magnetization at t = 0.
M_origin = [0; 0; 0];       % Origin of Magnetization vector.
scale = 1;

bloch = ode;
bloch.ODEFcn = @(t,M) [gamma*B0*M(2); -gamma*B0*M(1); 0*M(3)];
bloch.InitialTime = 0;
bloch.InitialValue = M_initial'; % Initial magnetization orientation. 
% Solve
solution = solve(bloch, time);
t = solution.Time;
M = solution.Solution;

% Plot x, y, and z components vs time.
figure
plot(t, M(1,:),'r','DisplayName',"Mx(t)");
hold on;
plot(t, M(2,:),'b','DisplayName',"My(t)");
hold on;
plot(t, M(3,:),'b','DisplayName',"Mz(t)");
legend();


% 3D Plot and appearance.
spinColors = [1 0 0]; 
showspins(M(:,1),scale,M_initial,spinColors);
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);

view(3)
title('3D View')
grid on
box on 
axis vis3d % Disable "stretch-to-fill" as view moves

% 3D Animation precession
for i=2:length(t)
    cla
    M_i = M(:,i);
    showspins(M_i,scale,M_origin, spinColors);
    pause(0.001);
    fprintf('time = %fs \n', t(i))
end


%% Free Precession in ROTATING Frame without Relaxation (Numerical Solver)

% Animation parameters.
Tend = 1;                   % Animation/Simulation duration (s).
frameRate = 60;             % # frames/second.
time = 0:1/frameRate:Tend;  % Animation/Simulation time points.
nFrames = Tend * frameRate; % Total frames.
gamma = 1;                  % Gyromagnetic ratio for H ([rad]/[s][T]).
B0 = 2*pi;                     % Tesla.
f = (gamma*B0)/(2*pi);      % Larmor frequency (Hz).
M_initial = [0; -0.5; 0.5]; % Magnetization at t = 0.
M_origin = [0; 0; 0];       % Origin of Magnetization vector.
scale = 1;

bloch = ode;
bloch.ODEFcn = @(t,M) [gamma*B0*M(2); -gamma*B0*M(1); 0*M(3)];
bloch.InitialTime = 0;
bloch.InitialValue = M_initial'; % Initial magnetization orientation. 
% Solve
solution = solve(bloch, time);
t = solution.Time;
M = solution.Solution;

% Plot x, y, and z components vs time.
figure
plot(t, M(1,:),'r','DisplayName',"Mx(t)");
hold on;
plot(t, M(2,:),'b','DisplayName',"My(t)");
hold on;
plot(t, M(3,:),'b','DisplayName',"Mz(t)");
legend();


% 3D Plot and appearance.
figure 
spinColors = [1 0 0]; 
showspins(M(:,1),scale,M_initial,spinColors);
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);

% set(gcf, 'Color', 'black');
% set(gca,'Color', 'black', ...
%         'XColor', 'black', ...
%         'YColor', 'black', ...
%         'ZColor', 'black')
view(3)
grid on
box on 
axis vis3d % Disable "stretch-to-fill" as view moves


% 3D Animation precession
viewDirection = zrot(-90)*M_initial; % View plane perpendicular to M_initial.
viewDirection = viewDirection(1:2);
initialView = acosd(dot(viewDirection,M_initial(1:2))/(norm(viewDirection)*norm(M_initial(1:2))));
view(initialView,30)
for i=2:length(t)
    cla
    M_i = M(:,i);
    showspins(M_i,scale,M_origin, spinColors);
    view(initialView, 20)
    initialView = initialView - f*(360/frameRate); % Rotate the view along with the rate of precession!
    pause(0.001);
    fprintf('time = %fs \n', t(i))
end


%% Free Precession in Laboratory Frame with Relaxation (Numerical Solver)

% Animation parameters.
Tend = 10;                   % Animation/Simulation duration (s).
frameRate = 60;             % # frames/second.
time = 0:1/frameRate:Tend;  % Animation/Simulation time points.
nFrames = Tend * frameRate; % Total frames.
gamma = 1;                  % Gyromagnetic ratio for H ([rad]/[s][T]).
B0 = 2*pi;                  % Tesla.
f = (gamma*B0)/(2*pi);      % Larmor frequency (Hz).
M_initial = [0; -1; 0]; % Magnetization at t = 0.
M_eq = [0; 0; 1];           % Equilibrium Magnetization.
M_origin = [0; 0; 0];       % Origin of Magnetization vector.
scale = 1;
T1 = 2;                     % T1 relaxation (s).
T2 = 1;                     % T2 relaxation (2).

bloch = ode;
bloch.ODEFcn = @(t,M) [gamma*B0*M(2) - M(1)/T2; -gamma*B0*M(1) - M(2)/T2; -(M(3) - M_eq(3))/T1];
bloch.InitialTime = 0;
bloch.InitialValue = M_initial'; % Initial magnetization orientation. 
% Solve
solution = solve(bloch, time);
t = solution.Time;
M = solution.Solution;

% Plot x, y, and z components vs time.
figure
plot(t, M(1,:),'r','DisplayName',"Mx(t)");
hold on;
plot(t, M(2,:),'b','DisplayName',"My(t)");
hold on;
plot(t, M(3,:),'b','DisplayName',"Mz(t)");
legend();


% 3D Plot and appearance.
figure
spinColors = [1 0 0]; 
showspins(M(:,1),scale,M_initial,spinColors);
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);

title('3D View')
grid on
box on 
axis vis3d % Disable "stretch-to-fill" as view moves

% 3D Animation precession
view(-35, 10)
for i=2:length(t)
    cla
    M_i = M(:,i);
    showspins(M_i,scale,M_origin, spinColors);
    pause(0.001);
    fprintf('time = %fs \n', t(i))
end


%% Forced Precession in Laboratory Frame without Relaxation (Numerical Solver)

% Animation parameters.
Tend = 10;                  % Animation/Simulation duration (s).
frameRate = 60;             % # frames/second.
time = 0:1/frameRate:Tend;  % Animation/Simulation time points.
nFrames = Tend * frameRate; % Total frames.
gamma = 1;                  % Gyromagnetic ratio for H ([rad]/[s][T]).
B0 = 2*pi;                  % Tesla.
w = gamma*B0;               % Larmor frequency(rad/s).
f = (gamma*B0)/(2*pi);      % Larmor frequency (Hz).
M_initial = [0; 0; 1];      % Magnetization at t = 0.
M_origin = [0; 0; 0];       % Origin of Magnetization vector.
scale = 3;                  % Scale the axis limits.

% B1 parameters: B1(t) = B1cos(wt)i - B1sin(wt)j (left handed rotation)
B1x = 1;  % Magnitude of x component of oscillating B1 
B1y = 1;  % Magnitude of y component of oscillating B1
B1z = 0;  % Magnitude of z component of oscillating B1
w1 = w;                     % Frequency of B1 pulse

% ODE Solver
bloch = ode;
bloch.ODEFcn = @(t,M) [gamma*M(2)*B0 + gamma*M(3)*B1y*sin(w1*t);
                      -gamma*M(1)*B0 + gamma*M(3)*B1x*cos(w1*t); 
                      -gamma*M(1)*B1y*sin(w1*t) - gamma*M(2)*B1x*cos(w1*t)];
bloch.InitialTime = 0;
bloch.InitialValue = M_initial'; % Initial magnetization orientation. 
% Solve
solution = solve(bloch, time);
t = solution.Time;
M = solution.Solution;

% Plot x, y, and z components vs time.
figure
componentAxes = axes;
plot(componentAxes, t, M(1,:),'r','DisplayName',"Mx(t)");
hold on;
plot(componentAxes, t, M(2,:),'b','DisplayName',"My(t)");
hold on;
plot(componentAxes, t, M(3,:),'b','DisplayName',"Mz(t)");
legend();

% 3D Plot and appearance.
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
spinColors = [1 0 0]; 
showspins(M(:,1),scale,M_origin,spinColors);
view(3)
grid on
box on 
axis vis3d % Disable "stretch-to-fill" as view moves
set(gcf, 'Color', 'black');
set(gca,'Color', 'black', ...
        'XColor', 'black', ...
        'YColor', 'black', ...
        'ZColor', 'black')

% 3D Animation precession
initialView = 90;
for i=2:length(t)
    cla
    M_i = M(:,i);
    showspins(M_i,scale,M_origin, spinColors);
    view(initialView, 20)
    initialView = initialView - f*(360/frameRate); % Rotate the view along with the rate of precession!
    
    % Plot axis lines
    line(2*xlim, [0,0], [0,0], 'LineWidth', 3, 'Color', 'b');
    line([0,0], 2*ylim, [0,0], 'LineWidth', 3, 'Color', 'b');
    line([0,0], [0,0], 2*zlim, 'LineWidth', 3, 'Color', 'b');
    pause(0.001);
    fprintf('time = %fs \n', t(i))
end
