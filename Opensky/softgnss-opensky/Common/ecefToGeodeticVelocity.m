function [velocity_geo] = ecefToGeodeticVelocity(v_ecef, X, Y, Z)

    % WGS84 ellipsoid constants
    a = 6378137; % Semi-major axis
    e2 = 6.69437999014e-3; % Square of eccentricity

    % Calculate longitude and latitude
    lon = atan2(Y, X);
    p = sqrt(X^2 + Y^2);
    theta = atan2(Z * a, p * (1 - e2));
    lat = atan2(Z + e2 * (a / sqrt(1 - e2) * sin(theta)^3), ...
                p - e2 * (a / sqrt(1 - e2)) * cos(theta)^3);

    % Calculate the transformation matrix
    N = a / sqrt(1 - e2 * sin(lat)^2); % Radius of curvature in the prime vertical

    % Jacobian matrix from ECEF to geodetic coordinates
    J = [-sin(lat) * cos(lon) / (N + Z), -sin(lon) / (N + Z), cos(lat) * cos(lon) / (N + Z);
          -sin(lat) * sin(lon) / (N + Z), cos(lon) / (N + Z), cos(lat) * sin(lon) / (N + Z);
          cos(lat) / (N * (1 - e2 * sin(lat)^2)), 0, -sin(lat) / (N * (1 - e2 * sin(lat)^2))];

    % Convert ECEF velocity to geodetic velocity
    velocity_geo = J * v_ecef; % v_geo = J * v_ecef

end