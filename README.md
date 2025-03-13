# AAE6102 Assignment 1

**Author:** CHENG Yuling  
**Student ID:** 24041314r  

## Overview
This repository contains the code and documentation for AAE Assignment 1. All tasks were completed using the softGNSS open-source code. The primary modifications were made in the following files:

- `postProcessing.m`
- `posNavigation.m`
- Additional plot code for drawing

All images in openskypng and urbanpng files are screenshots of the corresponding file tasks results.

## Task Details

### Task 1
Satellite signal acquisition under:
| Parameter              | Open-sky Dataset                        | Urban Dataset                  |
| ---------------------- | --------------------------------------- | ------------------------------ |
| Carrier frequency      | 1575.42 MHz                             | 1575.42 MHz                    |
| Intermediate frequency | 4.58 MHz                                | 0 MHz                          |
| Sampling frequency     | 58 MHz                                  | 26 MHz                         |
| Data format            | 8-bit I/Q samples                       | 8-bit I/Q samples              |
| Ground truth           | (22.328444770087565, 114.1713630049711) | (22.3198722, 114.209101777778) |
| Data length            | 90 seconds                              | 90 seconds                     |

#### Open Sky
| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |
|---------|-----|---------------|-----------|-------------|--------|
|    1    |  16 |  4.57976e+06 |    -240   |     31994   |   T    |
|    2    |  26 |  4.58192e+06 |    1917   |     57754   |   T    |
|    3    |  31 |  4.58107e+06 |    1066   |     18744   |   T    |
|    4    |  22 |  4.58157e+06 |    1571   |     55101   |   T    |
|    5    |  27 |  4.57678e+06 |   -3220   |      8814   |   T    |

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/aquicition.png alt="aquisition" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/navigation3.png alt="navigation" width=300 height=300>

#### Urban
| Channel | PRN  | Frequency   | Doppler | Code Offset | Status |
| ------- | ---- | ----------- | ------- | ----------- | ------ |
|     1   | 1   | 1.20258e+03 |    1203   |      3329   |     T  |
|     2   |   3 |  4.28963e+03 |    4290   |     25173   |     T  |
|     3   |  11 |  4.09126e+02 |     409   |      1155   |     T  |
|     4   |  18 |  -3.22342e+02 |    -322   |     10581   |     T  |

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/acquisition.png alt="acquisition" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/navigation3.png alt="navigation" width=300 height=300>

### Task 2
After acquisition, we can track
#### Open Sky
For example, for PRN 16, we have

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN16result.png alt="替代文本" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN16CNo.png alt="替代文本" width=400 height=300>

For 26, we have

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN26result.png alt="替代文本" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN26CNo.png alt="替代文本" width=400 height=300>

For correlation plots,

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN11cor.png alt="替代文本" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/PRN1cor.png alt="替代文本" width=400 height=300>

#### Urban
For PRN 1, we have

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/PRN1result.png alt="替代文本" width=400 height=300><img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/PRN1CNo.png alt="替代文本" width=400 height=300>

For correlation plots,

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/PRN1cor.png alt="替代文本" width=400 height=300>

In urban environments, signals from satellites may reflect off buildings and other structures, creating multiple signal paths and degrading the signal. The presence of multiple paths may result in multiple correlation peaks in signal processing, which can confuse signal processing algorithms, leading to errors in positioning and navigation.

### Task 3
Decode the navigation message and extract key parameters, we can have ephemeris data from eph.
#### Open Sky
| PRN      | 16          | 22          | 26          | 27          | 31          |
|----------|-------------|-------------|-------------|-------------|-------------|
| C_ic     | -1.01E-07   | -1.01E-07   | -2.05E-08   | 1.08E-07    | -1.14E-07   |
| omega_0  | -1.674261429| 1.272735322 | -1.812930701| -0.71747466 | -2.787272903|
| C_is     | 1.36E-07    | -9.31E-08   | 8.94E-08    | 1.15E-07    | -5.03E-08   |
| i_0      | 0.971603403 | 0.936454583 | 0.939912327 | 0.974727542 | 0.95588255  |
| C_rc     | 237.6875    | 266.34375   | 234.1875    | 230.34375   | 240.15625   |
| omega    | 0.679609497 | -0.887886686| 0.295685419 | 0.630881665 | 0.311626182 |
| omegaDot | -8.01E-09   | -8.67E-09   | -8.31E-09   | -8.02E-09   | -7.99E-09   |
| IODE_sf3 | 9           | 22          | 113         | 30          | 83          |
| iDot     | -4.89E-10   | -3.04E-11   | -4.18E-10   | -7.14E-13   | 3.21E-11    |
| idValid  | [2,0,3]     | [2,0,3]     | [2,0,3]     | [2,0,3]     | [2,0,3]     |
| weekNumber| 1155       | 1155        | 1155        | 1155        | 1155        |
| accuracy | 0           | 0           | 0           | 0           | 0           |
| health   | 0           | 0           | 0           | 0           | 0           |
| T_GD     | -1.02E-08   | -1.77E-08   | 6.98E-09    | 1.86E-09    | -1.30E-08   |
| IODC     | 234         | 218         | 15          | 4           | 228         |
| t_oc     | 396000      | 396000      | 396000      | 396000      | 396000      |
| a_f2     | 0           | 0           | 0           | 0           | 0           |
| a_f1     | -6.37E-12   | 9.21E-12    | 3.98E-12    | -5.00E-12   | -1.93E-12   |
| a_f0     | -0.000406925| -0.000489472| 0.00014479  | -0.000206121| -0.0001449  |
| IODE_sf2 | 9           | 22          | 113         | 30          | 83          |
| C_rs     | 23.34375    | -99.8125    | 21.25       | 70.4375     | 30.71875    |
| deltan   | 4.25E-09    | 5.28E-09    | 5.05E-09    | 4.03E-09    | 4.81E-09    |
| M_0      | 0.718116855 | -1.260965589| 1.735570934 | -0.173022281| 2.82452322  |
| C_uc     | 1.39E-06    | -5.16E-06   | 1.15E-06    | 3.73E-06    | 1.46E-06    |
| e        | 0.012296279 | 0.006713538 | 0.006253509 | 0.009574107 | 0.010271554 |
| C_us     | 7.69E-06    | 5.17E-06    | 7.04E-06    | 8.24E-06    | 7.23E-06    |
| sqrtA    | 5153.771322 | 5153.712273 | 5153.636459 | 5153.652021 | 5153.622389 |
| t_oe     | 396000      | 396000      | 396000      | 396000      | 396000      |
| TOW      | 390102      | 390102      | 390102      | 390102      | 390102      |

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/navigation.png alt="替代文本" width=600 height=500>

#### Urban
| PRN        | 1            | 3            | 11           | 18          |
| ---------- | ------------ | ------------ | ------------ | ----------- |
| C_ic       | -7.45E-08    | 1.12E-08     | -3.17E-07    | -2.53E-07   |
| omega_0    | -3.106035801 | -2.064178438 | 2.725770376  | 3.121821254 |
| C_is       | 1.60E-07     | 5.22E-08     | -1.32E-07    | 3.54E-08    |
| i_0        | 0.976127704  | 0.962858746  | 0.909806736  | 0.9546426   |
| C_rc       | 287.46875    | 160.3125     | 324.40625    | 280.15625   |
| omega      | 0.711497599  | 0.594974558  | 1.891492962  | 1.393015876 |
| omegaDot   | -8.17E-09    | -7.83E-09    | -9.30E-09    | -8.61E-09   |
| IODE_sf3   | 72           | 72           | 83           | 56          |
| iDot       | -1.81E-10    | 4.81E-10     | 1.29E-11     | -1.62E-10   |
| idValid    | [2,0,3]      | [2,0,3]      | [2,0,3]      | [2,0,3]     |
| weekNumber | 1032         | 1032         | 1032         | 1032        |
| accuracy   | 0            | 0            | 0            | 0           |
| health     | 0            | 0            | 0            | 0           |
| T_GD       | 5.59E-09     | 1.86E-09     | -1.26E-08    | -5.59E-09   |
| IODC       | 12           | 4            | 229          | 244         |
| t_oc       | 453600       | 453600       | 453600       | 453600      |
| a_f2       | 0            | 0            | 0            | 0           |
| a_f1       | -9.44E-12    | -1.14E-12    | 8.53E-12     | 3.18E-12    |
| a_f0       | -3.49E-05    | 0.000186326  | -0.000590093 | 5.99E-05    |
| IODE_sf2   | 72           | 72           | 83           | 56          |
| C_rs       | -120.71875   | -62.09375    | -67.125      | -113.875    |
| deltan     | 4.19E-09     | 4.45E-09     | 5.89E-09     | 4.72E-09    |
| M_0        | 0.517930888  | -0.430397464 | -0.198905418 | 0.259840989 |
| C_uc       | -6.33E-06    | -3.09E-06    | -3.60E-06    | -6.11E-06   |
| e          | 0.008923085  | 0.00222623   | 0.016643139  | 0.015419818 |
| C_us       | 5.30E-06     | 1.16E-05     | 1.51E-06     | 5.11E-06    |
| sqrtA      | 5153.655643  | 5153.777802  | 5153.706596  | 5153.699318 |
| t_oe       | 453600       | 453600       | 453600       | 453600      |
| TOW        | 449352       | 449352       | 449352       | 449352      |

<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/navigation.png alt="替代文本" width=600 height=400>

### Task 4
Consider a receiver that receives signals from multiple satellites and estimates its position using pseudoranges. The receiver's position and clock bias are estimated. In WLS, the weight of the error is considered and the following formula is used for parameter estimation: $\hat{\mathbf{p}} = (\mathbf{H}^T \mathbf{W} \mathbf{H})^{-1} \mathbf{H}^T \mathbf{W} \mathbf{r}$

The receiver's velocity can be calculated by taking the time derivative of the receiver's position estimate. In many cases, the receiver's velocity can be obtained by observing the change in the receiver's position at consecutive moments.
#### Open Sky
<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/Compare.png alt="替代文本" width=600 height=500>

#### Urban
<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/Compare.png alt="替代文本" width=600 height=500>

### Task 5

The Extended Kalman Filter (EKF) involves several key steps to estimate the state of a dynamic system:

#### Initialization:
   - Set the initial state estimate and covariance matrix.
   - Define the process noise and measurement noise covariance matrices.

#### Prediction Step:
   - Predict the current state based on the previous state estimate and control inputs.
   - Calculate the covariance of the state estimate, incorporating process noise.

#### Update Step:
   - Compute the measurement residual (the difference between the actual measurements and predicted measurements).
   - Calculate the Kalman gain to optimize the state update.
   - Update the state estimate by combining the predicted state and the measurement.
   - Update the covariance matrix to reflect the new uncertainty.

#### Iterative Processing:
   - Repeat the prediction and update steps to continuously process state estimates over time.

These steps enable the Extended Kalman Filter to effectively track the state of nonlinear dynamic systems while integrating measurements from sensors to provide real-time state estimations.
#### Open Sky
<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/openskypng/EKF.png alt="替代文本" width=400 height=300>

#### Urban
<img src=https://github.com/CHENG-Yuling/PolyU_AAE6102_Assignment1/blob/main/urbanpng/EKF.png alt="替代文本" width=400 height=300>


