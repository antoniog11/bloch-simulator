addpath '/Users/antonioglenn/Documents/code-and-projects/RAD229 MRI Signals and Sequences/Matlab'
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


%% Free/Forced Precession (Lab or Rotating Frame) with/without Relaxation (Numerical Solver)
close all

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
M_eq = [0; 0; 1];           % Equilibrium Magnetization.
scale = 1;
T1 = 3;                     % T1 relaxation (s).
T2 = 2;                     % T2 relaxation (2).
viewRotatingFrame = true;   % View precession in rotating reference frame
relaxationOn = false;       % Enable T1/T2 Relaxation?
viewElevation = 0;         % Elevation to view spins (degrees): view(~, viewElevation)
traceSpinTip = true;        % Enable tracing of tip of magnetization vector.

% B1 parameters: B1(t) = B1cos(wt)i - B1sin(wt)j (left handed rotation)
nTR = 6;
B1mag = 4;                  % Magnitude of B1
w1 = w;                     % Frequency of B1 pulse
B1phase = pi/2;             % Phase of B1 (radians)
B1phase_degrees = B1phase * (180/pi);
B1x = @(t) B1mag*cos(w1*t + B1phase); % x component of B1
B1y = @(t) B1mag*sin(w1*t + B1phase); % y component of B1
T_90 = (pi/2) / (gamma*B1mag);
T_180 = (pi/2) / (gamma*B1mag);
flipAngle = pi/8;

% ODE initialization
bloch = ode;
bloch.ODEFcn = @(t, M) blochODEFcn(t, M, gamma, B0, B1x, B1y, T1, T2, M_eq, relaxationOn, w1, B1mag, nTR, Tend, flipAngle);
bloch.InitialTime = 0;
bloch.InitialValue = M_initial; % Initial magnetization orientation. 
function dMdt = blochODEFcn(t, M, gamma, B0, B1x, B1y, T1, T2, M_eq, relaxationOn, w1, B1mag, nTR, Tend, flipAngle)
    % Calculate the time interval between excitations
    TR_interval = Tend / nTR;
    % Determine if t falls within an "on" phase of the B1 field
    B1Duration = flipAngle/(gamma*B1mag); % Duration neded for flip Angle of theta radians.
    % Modulo operation helps check the timing within each TR interval
    if mod(t, TR_interval) < B1Duration
        B1x_t = B1x(t);
        B1y_t = B1y(t);
    else
        B1x_t = 0;
        B1y_t = 0;
    end

    % Define the Bloch equations with or without relaxation
    if relaxationOn
        dMdt = [
            gamma * B0 * M(2) + gamma * M(3) * B1y_t - M(1) / T2;
            -gamma * B0 * M(1) + gamma * M(3) * B1x_t - M(2) / T2;
            -gamma * M(1) * B1y_t - gamma * M(2) * B1x_t - (M(3) - M_eq(3)) / T1
        ];
    else
        dMdt = [
            gamma * B0 * M(2) + gamma * M(3) * B1y_t;
            -gamma * B0 * M(1) + gamma * M(3) * B1x_t;
            -gamma * M(1) * B1y_t - gamma * M(2) * B1x_t
        ];
    end
end

% Solve
solution = solve(bloch, time);
t = solution.Time;
M = solution.Solution;

% Plot Magnetization components (x, y, and z) vs time.
figure
plot(t, M(1,:),'r','DisplayName',"Mx(t)");
hold on;
plot(t, M(2,:),'b','DisplayName',"My(t)");
hold on;
plot(t, M(3,:),'b','DisplayName',"Mz(t)");
legend();


% 3D Plot and appearance.
fig = figure(2);
fig.Position = [10 10 1000 1000];
%figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
spinColors = [1 0 0]; 
set(gcf, 'Color', 'black');
set(gca, 'XLim', [-2 2], ...                     % New Axes Limits
         'YLim', [-2 2], ...
         'ZLim', [-2 2]);
set(gca,'Color', 'black', ...
        'XColor', 'black', ...
        'YColor', 'black', ...
        'ZColor', 'black')
line(2*xlim, [0,0], [0,0], 'LineWidth', 3, 'Color', 'w');
line([0,0], 2*ylim, [0,0], 'LineWidth', 3, 'Color', 'w');
line([0,0], [0,0], 2*zlim, 'LineWidth', 3, 'Color', 'w');
view(3)
grid on
box on 
axis vis3d % Disable "stretch-to-fill" as view moves

if traceSpinTip
    spinTip = animatedline('linewidth', 2, 'Color', 'b');      % Animated Line tracing spin vector tip.
end

% 3D Animation precession
if viewRotatingFrame
    if B1mag ~= 0 % If we have a B1 field, align view along B1 axis
        initialView = B1phase_degrees - 90;
        view(initialView, viewElevation);
    else
        viewDirection = zrot(-90)*M_initial; % View plane perpendicular to M_initial.
        viewDirection = viewDirection(1:2);
        initialView = acosd(dot(viewDirection,M_initial(1:2))/(norm(viewDirection)*norm(M_initial(1:2))));
        view(initialView,viewElevation)
    end
end

spin = showspins(M(:,1),scale,M_initial,spinColors);
for i=1:length(t)
    figure(fig)
    delete(spin)
    spin = showspins(M(:,i), scale, M_origin, spinColors);
    if viewRotatingFrame
        view(initialView, viewElevation)
        initialView = initialView - f*(360/frameRate); % Rotate the view along with the rate of precession!
    end
      
    if traceSpinTip
        addpoints(spinTip, M(1,i), M(2,i), M(3,i));
    end

    pause(0.001);
    fprintf('time = %fs \n', t(i))
end