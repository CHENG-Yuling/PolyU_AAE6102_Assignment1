function [navSolutions, eph] = postNavigation(trackResults, settings)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis with help from Kristin Larson
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

%% Check is there enough data to obtain any navigation solution ===========
% It is necessary to have at least three subframes (number 1, 2 and 3) to
% find satellite coordinates. Then receiver position can be found too.
% The function requires all 5 subframes, because the tracking starts at
% arbitrary point. Therefore the first received subframes can be any three
% from the 5.
% One subframe length is 6 seconds, therefore we need at least 30 sec long
% record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
% when tracking has started in a middle of a subframe.

if (settings.msToProcess < 36000) 
    % Show the error message and exit
    disp('Record is to short. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Pre-allocate space =======================================================
% Starting positions of the first message in the input bit stream 
% trackResults.I_P in each channel. The position is PRN code count
% since start of tracking. Corresponding value will be set to inf 
% if no valid preambles were detected in the channel.
subFrameStart  = inf(1, settings.numberOfChannels);

% Time Of Week (TOW) of the first message(in seconds). Corresponding value
% will be set to inf if no valid preambles were detected in the channel.
TOW  = inf(1, settings.numberOfChannels);

%--- Make a list of channels excluding not tracking channels ---------------
activeChnList = find([trackResults.status] ~= '-');

%% Decode ephemerides =======================================================
for channelNr = activeChnList
    
    % Get PRN of current channel
    PRN = trackResults(channelNr).PRN;
    
    fprintf('Decoding NAV for PRN %02d -------------------- \n', PRN);
    %=== Decode ephemerides and TOW of the first sub-frame ==================
    [eph(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  NAVdecoding(trackResults(channelNr).I_P);  %#ok<AGROW>

    %--- Exclude satellite if it does not have the necessary nav data -----
    if (isempty(eph(PRN).IODC) || isempty(eph(PRN).IODE_sf2) || ...
        isempty(eph(PRN).IODE_sf3))

        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d !\n', PRN);
    else
        fprintf('    Three requisite messages for PRN %02d all decoded!\n', PRN);
    end
end

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Set measurement-time point and step  =====================================
% Find start and end of measurement point locations in IF signal stream with available
% measurements
sampleStart = zeros(1, settings.numberOfChannels);
sampleEnd = inf(1, settings.numberOfChannels);

for channelNr = activeChnList
    sampleStart(channelNr) = ...
          trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
end

% Second term is to make space to aviod index exceeds matrix dimensions, 
% thus a margin of 1 is added.
sampleStart = max(sampleStart) + 1;  
sampleEnd = min(sampleEnd) - 1;

%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod/1000);

%---  Number of measurment point from measurment start to end ------------- 
measNrSum = fix((sampleEnd-sampleStart)/measSampleStep);

%% Initialization =========================================================
% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.  
readyChnList = activeChnList;

% Set local time to inf for first calculation of receiver position. After
% first fix, localTime will be updated by measurement sample step.
localTime = inf;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################

fprintf('Positions are being computed. Please wait... \n');
for currMeasNr = 1:measNrSum

    fprintf('Fix: Processing %02d of %02d \n', currMeasNr,measNrSum);
    %% Initialization of current measurement ==============================          
    % Exclude satellites, that are belove elevation mask 
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
                              readyChnList);

    % Save list of satellites used for position calculation
    navSolutions.PRN(activeChnList, currMeasNr) = ...
                                        [trackResults(activeChnList).PRN]; 

    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions.az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
                                     
    % Signal transmitting time of each channel at measurement sample location
    navSolutions.transmitTime(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);
    navSolutions.satClkCorr(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);                                                                  
       
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
    
%% Find pseudoranges ======================================================
    % Raw pseudorange = (localTime - transmitTime) * light speed (in m)
    % All output are 1 by settings.numberOfChannels columme vecters.
    [navSolutions.rawP(:, currMeasNr),transmitTime,localTime]=  ...
                     calculatePseudoranges(trackResults,subFrameStart,TOW, ...
                     currMeasSample,localTime,activeChnList, settings);
    
    % Save transmitTime
    navSolutions.transmitTime(activeChnList, currMeasNr) = ...
                                        transmitTime(activeChnList);

%% Find satellites positions and clocks corrections =======================
    % Outputs are all colume vectors corresponding to activeChnList
    [satPositions, satClkCorr] = satpos(transmitTime(activeChnList), ...
                                 [trackResults(activeChnList).PRN], eph); 
                                    
    navSolutions.satPositions(:, :, currMeasNr) = satPositions; 
    % Save satClkCorr
    navSolutions.satClkCorr(activeChnList, currMeasNr) = satClkCorr;
    %%


%% Find receiver position =================================================
    % 3D receiver position can be found only if signals from more than 3
    % satellites are available  
    if size(activeChnList, 2) > 3

        %=== Calculate receiver position ==================================
        % Correct pseudorange for SV clock error
        clkCorrRawP = navSolutions.rawP(activeChnList, currMeasNr)' + ...
                                                   satClkCorr * settings.c;

        % Calculate receiver position
        [xyzdt,navSolutions.el(activeChnList, currMeasNr), ...
               navSolutions.az(activeChnList, currMeasNr), ...
               navSolutions.DOP(:, currMeasNr)] =...
                       leastSquarePos(satPositions, clkCorrRawP, settings);

        %=== Save results ===========================================================
        % Receiver position in ECEF
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);       
        
		% For first calculation of solution, clock error will be set 
        % to be zero
        if (currMeasNr == 1)
        navSolutions.dt(currMeasNr) = 0;  % in unit of (m)
        else
            navSolutions.dt(currMeasNr) = xyzdt(4);  
        end

        
                
		%=== Correct local time by clock error estimation =================
        localTime = localTime - xyzdt(4)/settings.c;       
        navSolutions.localTime(currMeasNr) = localTime;
        
        % Save current measurement sample location 
        navSolutions.currMeasSample(currMeasNr) = currMeasSample;

        % Update the satellites elevations vector
        satElev = navSolutions.el(:, currMeasNr)';

        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.correctedP(activeChnList, currMeasNr) = ...
                navSolutions.rawP(activeChnList, currMeasNr) + ...
                satClkCorr' * settings.c - xyzdt(4);
            
%% Coordinate conversion ==================================================

        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
         navSolutions.longitude(currMeasNr), ...
         navSolutions.height(currMeasNr)] = cart2geo(...
                                            navSolutions.X(currMeasNr), ...
                                            navSolutions.Y(currMeasNr), ...
                                            navSolutions.Z(currMeasNr), ...
                                            5);
        
        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
                                       navSolutions.longitude(currMeasNr));
        
        % Position in ENU
        [navSolutions.E(currMeasNr), ...
         navSolutions.N(currMeasNr), ...
         navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
                                                xyzdt(3), ...
                                                navSolutions.utmZone);

        %%
        %=== Find receiver position ===========================================
        % 假设 navSolutions 结构包含 X, Y, Z 向量用于存储位置
% 以及所需的观测数据和设计矩阵

if (currMeasNr > 1)
    % 提取之前和当前的观测数据，这里假设 obsData 是观测数据
    %obsData = navSolutions.observations(:, currMeasNr);
    
    % 提取之前的位置估计，用于计算新的位置估计
    previousPosition = [navSolutions.X(currMeasNr - 1), ...
                            navSolutions.Y(currMeasNr - 1), ...
                            navSolutions.Z(currMeasNr - 1)];
    currentPosition = [navSolutions.X(currMeasNr), ...
                           navSolutions.Y(currMeasNr), ...
                           navSolutions.Z(currMeasNr)];

    correctedP_ECEF = navSolutions.correctedP(1, currMeasNr);  % 提取第 currMeasNr 列的第 1 行（第1个卫星的伪距）

    % 构建设计矩阵 H
    H1 = [1, 0, 0;  % 对于 X
          0, 1, 0;  % 对于 Y
          0, 0, 1]; % 对于 Z

    % 计算伪距向量 b1，使用当前位置与伪距对应的点（假设伪距是实际距离）
    distance = norm(currentPosition - correctedP_ECEF);  % 当前位置与当前伪距的距离
    b1 = distance;  % 这里 b1 是一个标量值，表示实际测量的伪距

    % 构建 b1 为列向量（对于矩阵计算需要）
    b1 = b1(:);  % 确保 b1 是列向量


    % 定义权重矩阵（可根据实际情况调整权重）
    W1 = eye(size(b1, 1)); % 权重矩阵，如果有多个观测则调整尺寸

    % 使用加权最小二乘法求解位置
    try
        positionEstimate = (H1' * W1 * H1) \ (H1' * W1 * b1);  % 计算位置估计
    catch ME
        disp('Error in calculating position estimate:');
        disp(ME.message);
        return;  % 如果计算过程中有误，退出
    end

    % 保存接收器位置
    navSolutions.X1(currMeasNr) = positionEstimate(1);
    navSolutions.Y1(currMeasNr) = positionEstimate(2);
    navSolutions.Z1(currMeasNr) = positionEstimate(3);
    [navSolutions.latitude1(currMeasNr), ...
         navSolutions.longitude1(currMeasNr), ...
         navSolutions.height1(currMeasNr)] = cart2geo(...
                                            navSolutions.X1(currMeasNr), ...
                                            navSolutions.Y1(currMeasNr), ...
                                            navSolutions.Z1(currMeasNr), ...
                                            5);

else
    % 第一个测量点的位置信息设置为 NaN
    navSolutions.X1(currMeasNr) = NaN;
    navSolutions.Y1(currMeasNr) = NaN;
    navSolutions.Z1(currMeasNr) = NaN;
end
        %%
   if (currMeasNr > 1)
        % 计算位置变化的最小二乘法
        previousPosition = [navSolutions.X(currMeasNr - 1), ...
                            navSolutions.Y(currMeasNr - 1), ...
                            navSolutions.Z(currMeasNr - 1)];
        currentPosition = [navSolutions.X(currMeasNr), ...
                           navSolutions.Y(currMeasNr), ...
                           navSolutions.Z(currMeasNr)];

        timeInterval = settings.navSolPeriod / 1000; % 单位：秒

        % 构建设计矩阵 H 和伪距向量 b
        H = [1, 0, 0; 0, 1, 0; 0, 0, 1]; % 设计矩阵
        b = (currentPosition - previousPosition) / timeInterval; % 速度变化

        % 使用单位权重矩阵（可根据实际情况调整权重）
        W = eye(3); % 权重矩阵，3x3单位矩阵

        % 使用加权最小二乘法求解速度
        velocity = (H' * W * H) \ (H' * W * b'); % b' 转置为列向量

        % 保存接收器速度
        navSolutions.velocity(:, currMeasNr) = velocity;
        velocity1 = navSolutions.velocity(:, currMeasNr); % 从 your navSolutions 获取速度
        X2 = navSolutions.X(currMeasNr); % ECEF X 坐标
        Y2 = navSolutions.Y(currMeasNr); % ECEF Y 坐标
        Z2 = navSolutions.Z(currMeasNr); % ECEF Z 坐标

        % 转换到大地坐标系
       velocity_geo= ecefToGeodeticVelocity(velocity1, X2, Y2, Z2);
       navSolutions.velocity_geo(:, currMeasNr) =velocity_geo;


    else
        navSolutions.velocity(:, currMeasNr) = NaN(3, 1); % 第一个测量点的速度设置为 NaN
    end



        
    else
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
                       ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;

        navSolutions.az(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
        navSolutions.el(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));

        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.

    end % if size(activeChnList, 2) > 3

    %=== Update local time by measurement  step  ====================================
    localTime = localTime + measSampleStep/settings.samplingFreq ;

end %for currMeasNr...
