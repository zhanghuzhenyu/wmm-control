# Adaptive Backstepping Tracking Control of Wheeled Mobile Manipulator (WMM) with Output Constraints

In this study, an adaptive tracking control scheme with output constraints is proposed for wheeled mobile manipulators. Firstly, a backstepping controller utilizing a time-varying barrier Lyapunov function is proposed to limit the tracking error within prescribed bounds. Further, a radial basis function neural network is employed to compensate for dynamic uncertainties and external disturbances.

## Code structure
- [WMM system](# WMM system)
- [Adaptive controller](# Adaptive controller)
- [Simulink](# Simulink)
- [Plot results](# Plot results)

## WMM system
The dynamic matrices of the WMM system are located in `sys.m`, which is a s-function.

## Adaptive controller
The adaptive controller of the WMM system is located in `controller.m`, which is a s-function.

## Simulink
We used Simulink to connect two parts `sys.m` and `controller.m` in `chap4_2sim.mdl`.

## Plot results
You can use `plot_result.m` to plot the tracking results after running simulation.
