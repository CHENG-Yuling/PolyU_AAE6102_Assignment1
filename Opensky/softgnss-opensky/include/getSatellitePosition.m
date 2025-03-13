function satellitePositions = getSatellitePosition(PRN, eph)
    % Constants
    mu = 3.986005e14; % Earth's gravitational parameter (m^3/s^2)
    omegaE = 7.2921151467e-5; % Earth's rotation rate (rad/s)

    % 从 eph 中提取参数
    % 根据 PRN 获取对应的卫星参数
    % 假设 eph 是一个结构数组，PRN 对应的索引为 PRN-1
    A = eph(PRN).sqrtA^2; % 半长轴
    e = eph(PRN).e; % 偏心率
    i0 = eph(PRN).i_0; % 升交点初始角
    omega = eph(PRN).omega; % 近地点角
    M0 = eph(PRN).M_0; % 平近点角
    Omega = eph(PRN).omega_0; % 升交点经度
    TOW = eph(PRN).TOW; % 周内时间

    % 计算卫星的平均运动
    n = sqrt(mu / A^3); % 平均运动 (rad/s)

    % 计算当前的平近点角
    M = M0 + n * TOW; % 当前的平近点角

    % 解决平近点方程
    E = M; % 初始猜测
    for j = 1:10
        E = M + e * sin(E); % 迭代求解
    end

    % 计算真近点角
    nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

    % 计算卫星的轨道位置
    r = A * (1 - e * cos(E)); % 距离
    x = r * (cos(Omega) * cos(omega + nu) - sin(Omega) * sin(omega + nu) * cos(i0));
    y = r * (sin(Omega) * cos(omega + nu) + cos(Omega) * sin(omega + nu) * cos(i0));
    z = r * (sin(i0) * sin(omega + nu));

    % 组合位置
    satellitePositions = [x, y, z]; % 卫星的位置 (x, y, z)
end