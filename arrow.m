% ~ arrow ~ >
% Silas Henderson (2025). Arrow 3d (https://www.mathworks.com/matlabcentral/fileexchange/71994-arrow-3d), MATLAB Central File Exchange. Retrieved February 5, 2025.
% Modified arrow function to accept list of colors.
function a = arrow(x0, y0, z0, x1, y1, z1, colors)
    X = [x0, x1, nan]; Y = [y0, y1, nan]; Z = [z0, z1, nan];
    u = [x1 - x0; y1 - y0; z1 - z0]; uN = norm(u);
    for i = 1:200
        rad = cross(u, ones(3, 1) - 2*rand(3, 1)); 
        rad = rad/norm(rad)*uN/20;   
        X = [X, x1, x1 - (x1 - x0)/5 + rad(1), nan];
        Y = [Y, y1, y1 - (y1 - y0)/5 + rad(2), nan];
        Z = [Z, z1, z1 - (z1 - z0)/5 + rad(3), nan];
    end
    a = line(X, Y, Z, 'linewidth', 3, 'color', colors); end
