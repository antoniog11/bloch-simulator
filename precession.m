%% Dependencies
% arrow.m on Matlab File Exchange
% showspins.m from Brian Hargreaves and Daniel Ennis RAD229 Course https://github.com/mribri999/MRSignalsSeqs 
addpath '/Users/antonioglenn/Documents/code-and-projects/RAD229 MRI Signals and Sequences/Matlab'

% Free/Forced Precession (Lab or Rotating Frame) with/without Relaxation (Numerical Solver)
close all

%% USER CONTROLLED PARMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters controlled by you to changed the simulation and animation.
% Spin parameters
B0 = 2*pi;                  % B0 field strenght [Tesla]
gamma = 1;                  % Gyromagnetic ratio for Hydrogen [rad]/[s][T]
M_initial = [0; 0; 1];      % Initial Magnetization vector at t = 0 [x; y; z]
Tend = 10;                  % Simulation duration [s]
T1 = 3;                     % T1 relaxation [s]
T2 = 2;                     % T2 relaxation [s]
% B1 pulse parameters
nTR = 1;                    % Number of B1 excitations
B1mag = 4;                  % Magnitude of B1
B1phase = pi/2;             % Phase of B1 (radians)
w1 = gamma*B0;              % Frequency of B1 pulse [radians/s]
flipAngle = 20;             % B1 flip angle [radians]
% Animation parameters
rotatingFrameOn = true;     % View precession in rotating reference frame?
relaxationOn = false;       % Enable T1/T2 Relaxation?
traceSpinTipOn = false;      % Enable tracing of tip of magnetization vector?
viewElevation = 0;          % Elevation to view spins [degrees] -> view(~, viewElevation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation 
% Animation parameters.
frameRate = 60;             % # frames/second.
time = 0:1/frameRate:Tend;  % Animation/Simulation time points.
nFrames = Tend * frameRate; % Total frames.
w = gamma*B0;               % Larmor frequency(rad/s).
f = (gamma*B0)/(2*pi);      % Larmor frequency (Hz).
M_origin = [0; 0; 0];       % Origin of Magnetization vector.
M_eq = [0; 0; 1];           % Equilibrium Magnetization.
scale = 1;
% B1 parameters: B1(t) = B1cos(wt)i - B1sin(wt)j (left handed rotation)
B1x = @(t) B1mag*cos(w1*t + B1phase); % x component of B1
B1y = @(t) B1mag*sin(w1*t + B1phase); % y component of B1
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
fig1 = figure(1);
ax1 = subplot(2,1,1);
hold(ax1, "on")
plot(ax1, t, M(1,:),'r','DisplayName',"Mx(t)");
plot(ax1, t, M(2,:),'b','DisplayName',"My(t)");
title("Transverse Magnetization Components")
xlabel("time (s)");
ylabel("Magnitude")
legend()
ax2 = subplot(2,1,2);
hold(ax2, "on")
plot(ax2, t, M(3,:),'b','DisplayName',"Mz(t)");
% Example of adding xlines to each subplot independently
line1 = xline(ax1, 0, 'black', 'LineWidth', 2); % Add xline to the top subplot
line2 = xline(ax2, 0, 'black', 'LineWidth', 2); % Add xline to the bottom subplot
title("Longitudinal Magnetization Components")
xlabel("time (s)");
ylabel("Magnitude")
legend();
% 3D Plot and appearance.
fig2 = figure(2);
fig2.Position = [10 10 1000 1000];
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
if traceSpinTipOn
    spinTip = animatedline('linewidth', 2, 'Color', 'b');      % Animated Line tracing spin vector tip.
end
% Rotating frame vs Laboratory frame setup
if rotatingFrameOn
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
% 3D Plotting
spin = showspins(M(:,1),scale,M_initial,spinColors);
for i=1:length(t)
    figure(fig2)
    delete(spin)
    spin = showspins(M(:,i), scale, M_origin, spinColors);
    if rotatingFrameOn
        view(initialView, viewElevation)
        initialView = initialView - f*(360/frameRate); % Rotate the view along with the rate of precession!
    end
      
    if traceSpinTipOn
        addpoints(spinTip, M(1,i), M(2,i), M(3,i));
    end

    pause(0.001);
    % TO DO: Parallelize the tracing of xlines on fig1.
    % figure(fig1)
    % line1.Value = t(i);
    % line2.Value = t(i);
    fprintf('time = %fs \n', t(i))
end