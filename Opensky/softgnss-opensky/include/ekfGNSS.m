function [estimatedStatesHistory, P] = ekfGNSS(measurements, initialState, Q, R)
    % ekfGNSS - 扩展卡尔曼滤波器用于 GNSS 状态估计
    %
    % Inputs:
    %   measurements - 包含伪距和多普勒信息的结构数组 (当前不使用)
    %   initialState - 状态矩阵 [状态1, 状态2, ...] 每列一个状态
    %   Q - 过程噪声协方差矩阵
    %   R - 测量噪声协方差矩阵
    %
    % Outputs:
    %   estimatedStatesHistory - 估计后状态的历史 [位置; 速度] 6 x numMeasurements
    %   P - 协方差矩阵

    % 状态向量维度
    [n, numMeasurements] = size(initialState);  % n 应为 6
    if n ~= 6
        error('initialState must have 6 rows.');
    end
    
    % 初始化状态和协方差矩阵
    estimatedStates = initialState(:, 1);  % 取第一列作为初始状态
    P = eye(n);  % 初始化协方差矩阵为单位矩阵

    % 初始化历史状态数组
    estimatedStatesHistory = zeros(n, numMeasurements);  % 6 x numMeasurements

    % 状态转移矩阵 (假设状态在短时间内保持不变)
    F = eye(n);  % 状态转移矩阵
    
    % 固定的伪距和多普勒值
    fixedPseudorange = 2.015069590964406e+07;  % 例如，设置固定的伪距值
    fixedDoppler = -240;         % 例如，设置固定的多普勒值

    % 循环处理每个测量
    for k = 1:numMeasurements
        % 预测的状态
        estimatedStatesPred = F * estimatedStates;  % 预测的状态
        P = F * P * F' + Q;  % 更新协方差

        % 使用固定的伪距和多普勒值
        pseudorange = fixedPseudorange;
        doppler = fixedDoppler;

        % 计算测量预测值
        predictedPseudorange = norm(estimatedStatesPred(1:3));  % 位置为状态向量的前3个元素
        predictedDoppler = estimatedStatesPred(4:6);  % 速度为状态向量的后3个元素

        % 雅可比矩阵 H
        H = [predictedPseudorange, 0, 0, 0, 0, 0;  % 伪距的雅可比矩阵
             0, 0, 0, predictedDoppler(1), predictedDoppler(2), predictedDoppler(3)]; % 多普勒的雅可比矩阵

        % 残差
        z = [pseudorange; doppler];  % 测量向量
        z_hat = [predictedPseudorange; norm(predictedDoppler)];
        y = z - z_hat;  % 残差

        % 协方差更新
        S = H * P * H' + R;  % 残差协方差
        K = P * H' / S;  % 卡尔曼增益



        % 更新状态
        estimatedStates = estimatedStatesPred + K * y;  % 更新状态
        P = (eye(n) - K * H) * P;  % 更新协方差

        % 存储当前的状态估计
        estimatedStatesHistory(:, k) = estimatedStates;  % 存入历史记录

        % 如果 k < numMeasurements, 取下一个状态作为初始状态
        if k < numMeasurements
            estimatedStates = initialState(:, k + 1);  % 更新为下一个时间步骤的初始状态
        end
    end
end