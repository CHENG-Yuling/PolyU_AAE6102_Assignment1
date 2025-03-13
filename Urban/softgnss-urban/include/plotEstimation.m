function plotEstimation(navSolutions, settings)
    % Function to plot estimated position vs ground truth.
    % Inputs:
    %  - navSolutions: structure containing estimated positions and velocities
    %  - settings: structure containing true position and time vector

    % Extract estimated position from navSolutions (longitude, latitude, height)
    estimatedLongitude = navSolutions.longitude;  % 1x179
    estimatedLatitude = navSolutions.latitude;    % 1x179
    estimatedHeight = navSolutions.height;        % 1x179

    % Combine longitude, latitude, and height into a 179x3 matrix
    estimatedPosition = [estimatedLongitude', estimatedLatitude', estimatedHeight'];

    % Extract estimated velocity from navSolutions
    estimatedVelocity = navSolutions.velocity_geo;    % 3x179

    % Add time vector based on navSolPeriod and number of velocity samples
    numSamples = 179;  % 假设您有178个样本
    time = (0:numSamples-1) * (settings.navSolPeriod / 1000);  % 转换为秒

    % Validate input dimensions
    if size(estimatedPosition, 2) ~= 3 || size(estimatedVelocity, 1) ~= 3
        error('Position must be Nx3 matrix and velocity must be 3xN matrix.');
    end

    % Extract ground truth position from settings
    groundTruthPosition = [
        settings.truePosition.E , ...
        settings.truePosition.N , ...
        settings.truePosition.U 
    ];

    % Validate time vector length
    if length(time) ~= size(estimatedVelocity, 2)
        error('Time vector length must match the number of velocity samples.');
    end

    % Plot User Position vs Ground Truth
    figure;

    % Position plot
    subplot(2, 1, 1);
    plot3(estimatedPosition(:, 1), estimatedPosition(:, 2), estimatedPosition(:, 3), 'b+', 'DisplayName', 'Estimated Position');
    hold on;
    plot3(groundTruthPosition(:, 1), groundTruthPosition(:, 2), groundTruthPosition(:, 3), 'r+', 'LineWidth', 1.5, 'MarkerSize', 10, 'DisplayName', 'Ground Truth Position');
    grid on;
    xlabel('Longitude');
    ylabel('Latitude');
    zlabel('Height');
    title('Estimated Position vs Ground Truth');
    legend show;

    % Velocity plot
    subplot(2, 1, 2);
    plot(time, estimatedVelocity(1, :), 'b-', 'DisplayName', 'Estimated Velocity X');
    hold on;
    plot(time, estimatedVelocity(2, :), 'g-', 'DisplayName', 'Estimated Velocity Y');
    plot(time, estimatedVelocity(3, :), 'r-', 'DisplayName', 'Estimated Velocity Z');
    grid on;
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Estimated Velocity');
    legend show;

    % Adjust layout
    sgtitle('User Position and Velocity Comparison');
end