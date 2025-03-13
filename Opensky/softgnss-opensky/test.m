A = [0.569754283638669	0.0421158879058621	0.820735224211370	-1
0.592581550125153	0.554537892705594	0.584238677258628	-1
0.849968313695081	0.439366868073725	-0.290707105096223	-1
0.455410057942118	-0.472904718799163	-0.754296232300441	-1
0.769438978511893	-0.237166903071947	0.593056083717077	-1];

% 处理 NaN 值
A(isnan(A)) = 0;

P = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];

b = [-5707586.21759784
-4116223.03188849
-3315205.39577162
-3979173.04503699
-4523535.04809618];

% 添加正则化
lambda = 1e-5;  % 选择一个小的正数
x = (A' * P * A + lambda * eye(size(A, 2))) \ (A' * P * b);

%%
dopplerMeasurements = zeros(settings.numberOfChannels, 1);
satellitePositions = zeros(settings.numberOfChannels, 3); % 3D 坐标

% 提取伪距和卫星位置
for channelNr = 1:settings.numberOfChannels
    if trackResults(channelNr).PRN ~= 0
        % 计算Doppler测量
        dopplerMeasurements(channelNr) = channel(channelNr).acquiredFreq - settings.IF; 
        
        % 提取卫星位置
        satellitePositions(channelNr, :) = getSatellitePosition(trackResults(channelNr).PRN, eph);
    end
end

% 有效卫星数量
validSatellites = sum(any(satellitePositions, 2)); % 计算有效卫星数量
A = zeros(validSatellites, 4);
b = zeros(validSatellites, 1);

validIndex = 1; % 用于跟踪有效卫星的索引

% 获取接收机位置
%receiverPosition = [navSolutions.E, navSolutions.N, navSolutions.U]; % 使用 navSolutions 中的接收机位置

    for i = 1:size(satellitePositions, 1)
        if all(satellitePositions(i, :) ~= 0)
            satPos = satellitePositions(i, :);
            rho = navSolutions.correctedP(i, end); % 使用修正后的伪距
            if rho > 0
                receiverPosition = [navSolutions.X(end), navSolutions.Y(end), navSolutions.Z(end)];
                d = norm(satPos - receiverPosition);
                if d > 0
                    A(validIndex, :) = [(satPos(1) - receiverPosition(1)) / d, ...
                                        (satPos(2) - receiverPosition(2)) / d, ...
                                        (satPos(3) - receiverPosition(3)) / d, ...
                                        -1];
                
                    b(validIndex) = rho - d;
                    validIndex = validIndex + 1;
                else
                    disp(['Distance d is not valid for satellite ', num2str(i)]);
                end
            else
                disp(['Invalid rho for satellite ', num2str(i)]);
            end
        else
            disp(['Invalid satellite position for channel ', num2str(i)]);
        end
    end

% 清理 A 矩阵中的 NaN 值
A(isnan(A)) = 0;

% 加权矩阵 P（可以根据实际情况设置权重）
P = eye(validSatellites);
lamda = 1e-5; % 选择一个小的正数以提高数值稳定性

% 进行加权最小二乘法计算
x = (A' * P * A + lamda * eye(size(A, 2))) \ (A' * P * b);

% 提取估计的位置和速度
estimatedPosition = x(1:3); % 估计的位置
estimatedVelocity = x(4); % 如果需要，根据具体实现进行计算
%%
%% Calculate navigation solutions =========================================

% 准备 EKF 所需的参数
% 提取导航解中的估计位置（经度、纬度、高度）
estimatedLongitude = navSolutions.longitude;  % 1x178
estimatedLatitude = navSolutions.latitude;    % 1x178
estimatedHeight = navSolutions.height;        % 1x178

% 将经度、纬度和高度组合成 178x3 矩阵
estimatedPosition = [estimatedLongitude', estimatedLatitude', estimatedHeight']; % 178x3 矩阵

% 提取导航解中的估计速度
estimatedVelocity = navSolutions.velocity;    % 3x178
activeChnList = find([trackResults.status] ~= '-');

% 假设初始状态为导航解
% 这里我们需要将初始状态调整为合适的格式，以便 EKF 使用。
initialState = [estimatedPosition(1, :)'; estimatedVelocity(:, 7)]; % 取第一个时间步的状态
Q = diag([0.1, 0.1, 0.1, 0.01, 0.01, 0.01]); % 过程噪声协方差
R = diag([1, 0.1]); % 测量噪声协方差

% 假设您有一个总的测量数量
totalMeasurements = 178; % 根据实际情况调整

% 在处理每一批测量时，设置 currMeasNr
for currMeasNr = 1:totalMeasurements
    % 创建测量数据结构
    measurements = struct();
    for k = 1:length(activeChnList)
          channelNr = activeChnList(k);
       
            measurements(k).pseudorange = navSolutions.correctedP(channelNr, currMeasNr); % 提取对应通道的伪距
            
            % 计算 Doppler 值
            measurements(k).Doppler = channel(channelNr).acquiredFreq - settings.IF;
          
    end
    % 调用 EKF
    [estimatedStates, P] = ekfGNSS(measurements, initialState, Q, R);
    
    % 可以在此处处理或保存每个 currMeasNr 的结果
end
%%
% Example inputs
receiver_pos = [1000; 2000; 3000];  % Receiver position in ECEF (meters)
receiver_vel = [10; 20; 30];        % Receiver velocity in ECEF (m/s)
sat_pos = [20000; 10000; 15000];    % Satellite position in ECEF (meters)
sat_vel = [100; 200; 300];          % Satellite velocity in ECEF (m/s)
pseudorange = 22000;                % Measured pseudorange (meters)
doppler = 150;                      % Measured Doppler (Hz)
clock_bias = 0.1;                   % Receiver clock bias (meters)
clock_drift = 0.01;                 % Receiver clock drift (m/s)

% Calculate residuals
[delta_rho, delta_d] = calculate_residuals(receiver_pos, receiver_vel, sat_pos, sat_vel, pseudorange, doppler, clock_bias, clock_drift);

% EKF update (assuming H and R are already defined)
[x_updated, P_updated] = ekf_update(x_pred, P_pred, delta_rho, delta_d, H, R);

%% Calculate navigation solutions =========================================
disp('   Calculating navigation solutions...');

[navSolutions, eph] = postNavigation(trackResults, settings);

% Prepare for EKF processing
estimatedLongitude = navSolutions.longitude;  % 1x178
estimatedLatitude = navSolutions.latitude;    % 1x178
estimatedHeight = navSolutions.height;        % 1x178

% Combine longitude, latitude, and height into a 178x3 matrix
estimatedPosition = [estimatedLongitude; estimatedLatitude; estimatedHeight]; % 3x178 matrix


activeChnList = find([trackResults.status] ~= '-');

% Extract estimated velocity
estimatedVelocity = navSolutions.velocity;    % 3x178
estimatedVelocity(isnan(estimatedVelocity)) = 0;
initialState = [estimatedPosition; estimatedVelocity]; % Take the state from the first timestep

% Process noise covariance
Q = diag([0.1, 0.1, 0.1, 0.01, 0.01, 0.01]); % Q matrix
R = diag([1, 0.1]); % R matrix

% Number of measurements
totalMeasurements = size(estimatedPosition', 1); % 178 (or however many measurements exist)

% Initialize arrays to store estimated states
%estimatedStatesHistory = zeros(6, totalMeasurements); % 6 states (x, y, z, vx, vy, vz)

for currMeasNr = 1:totalMeasurements
    % Create measurement structure
    measurements = struct();
    
    % For each active channel, fill in pseudorange and Doppler measurements
    for k = 1:min(length(activeChnList), 5)  % Limit to 5 channels
        channelNr = activeChnList(k);
       
        % Extract pseudorange
        measurements.pseudorange = navSolutions.correctedP(1, currMeasNr);
        
        % Calculate and assign Doppler
        observedFreq = channel(1).acquiredFreq;  % Extract the received frequency
        IF = settings.IF;  % Intermediate frequency
        measurements.Doppler = observedFreq - IF;  % Calculate Doppler

 
    end

    % Call EKF (you would need to define the ekfGNSS function)
    [estimatedStates, P] = ekfGNSS(measurements, initialState, Q, R);
 
end
