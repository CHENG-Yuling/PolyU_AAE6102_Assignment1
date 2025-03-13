function plotEKF(estimatedStates)
% 假设 estimatedStates 的维度是 6 x 178
% 提取位置和速度
estimatedLongitude = estimatedStates(1, :); % 经度
estimatedLatitude = estimatedStates(2, :);  % 纬度
estimatedHeight = estimatedStates(3, :);    % 高度
estimatedVelocityX = estimatedStates(4, :); % 速度 X
estimatedVelocityY = estimatedStates(5, :); % 速度 Y
estimatedVelocityZ = estimatedStates(6, :); % 速度 Z

% 绘制位置和速度
figure;

% 绘制位置
subplot(2, 1, 1);
plot(estimatedLongitude, 'r', 'DisplayName', 'Longitude');
hold on;
plot(estimatedLatitude, 'g', 'DisplayName', 'Latitude');
plot(estimatedHeight, 'b', 'DisplayName', 'Height');
title('Estimated Position');
xlabel('Time Step');
ylabel('Position');
legend show;
grid on;

% 绘制速度
subplot(2, 1, 2);
plot(estimatedVelocityX, 'r', 'DisplayName', 'Velocity X');
hold on;
plot(estimatedVelocityY, 'g', 'DisplayName', 'Velocity Y');
plot(estimatedVelocityZ, 'b', 'DisplayName', 'Velocity Z');
title('Estimated Velocity');
xlabel('Time Step');
ylabel('Velocity (m/s)');
legend show;
grid on;

% 如果您有 plotEKF 函数，调用它

% 显示图形
sgtitle('EKF Estimated States');