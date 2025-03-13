function estimatedStates = ekf_positioning(pseudoranges, dopplerMeasurements, satellitePositions, initialState, processNoise, measurementNoise)
    % 扩展卡尔曼滤波器 (EKF) 实现
    % 输入:
    %   pseudoranges - 伪距测量
    %   dopplerMeasurements - 多普勒测量
    %   satellitePositions - 卫星位置
    %   initialState - 初始状态 [X0; Y0; Z0; Vx0; Vy0; Vz0]
    %   processNoise - 过程噪声协方差
    %   measurementNoise - 观测噪声协方差
    % 输出:
    %   estimatedStates - 估计的状态矩阵

    dt = 1; % 时间间隔，单位：秒
    numSatellites = length(pseudoranges);
    
    % 初始化状态向量
    state = initialState; % [X0; Y0; Z0; Vx0; Vy0; Vz0]
    
    % 初始化协方差矩阵
    O = eye(6); % 初始协方差矩阵
    Q = processNoise; % 过程噪声协方差
    R = measurementNoise; % 观测噪声协方差

    % 存储估计结果
    estimatedStates = zeros(6, length(pseudoranges));

    for k = 1:length(pseudoranges)
        % 预测步骤
        F = [1, 0, 0, dt, 0, 0;   % 状态转移矩阵
             0, 1, 0, 0, dt, 0;
             0, 0, 1, 0, 0, dt;
             0, 0, 0, 1, 0, 0;
             0, 0, 0, 0, 1, 0;
             0, 0, 0, 0, 0, 1];
        
        state = F * state; % 状态预测
        O = F * O * F' + Q; % 协方差预测

        % 观测步骤
        Z = [pseudoranges(k); dopplerMeasurements(k)]; % 当前伪距和多普勒测量
        H = zeros(numSatellites, 6); % 观测矩阵
        for i = 1:numSatellites
            satPos = satellitePositions(i, :);
            rho = norm(satPos - state(1:3)'); % 计算卫星到用户的距离
            H(i, :) = [(satPos(1) - state(1)) / rho, ...
                        (satPos(2) - state(2)) / rho, ...
                        (satPos(3) - state(3)) / rho, ...
                        0, 0, 0]; % 只更新位置部分
        end

        % 计算卡尔曼增益
        K = O * H' / (H * O * H' + R);
        
        % 更新步骤
        predictedMeasurement = H * state; % 预测测量值
        state = state + K * (Z - predictedMeasurement); % 更新状态
        O = (eye(6) - K * H) * O; % 更新协方差

        % 存储估计值
        estimatedStates(:, k) = state;
    end
end