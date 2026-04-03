# GBS3_trimmed1

This set corresponds to:

- `aircraft/GBS3_trimmed1/GBS3_trimmed1.xml`
- `aircraft/GBS3_trimmed1/GBS3_trimmed1_init.xml`
- `scripts/GBS3_trimmed1_script.xml`
- `engine/Zenoah_G-26A_GBS3_trimmed1.xml`

## Changes with respect to the original

### Inertia and mass

- `Ixz`: `0.0951 -> 0.08338`
- `CG.x`: `14.4881 in -> 19.9000 in`

### Aerodynamics

Only these two drag terms change:

- `CDo`: `0.1000 -> 0.053864`
- `Drag_induced`: `0.0877 -> 0.014951`

### Propulsion and FCS

- Engine: `Zenoah_G-26A -> Zenoah_G-26A_GBS3_trimmed1`
- Power: `1789.68 W -> 121.75 W`
- The throttle channel 

### Init

`GBS3_trimmed1_init.xml`:

- `vt = 20 kt`
- `altitude = 1000 ft`
- `phi = -0.11 deg`
- `theta = 7.81 deg`
- `alpha = 7.81 deg`

### Script

`GBS3_trimmed1_script.xml`:

- total duration: `5 s`
- `dt = 0.01`
- keeps the seeded controls until `t < 1.0 s`
- performs an early stabilizing trim at `t = 0.05 s`
- performs the full trim at `t = 1.0 s`

Control:

- `aileron = -0.010`
- `rudder = 0.000`
- `pitch-trim = -0.04`
- `throttle = 1.000000`

## Summary

`GBS3_trimmed1` keeps the original geometry and almost all the original aerodynamics. The changes with respect to the original model are:

- `Ixz`
- `CG`
- `CDo`
- `Drag_induced`
- engine power
- throttle channel 
- trim-specific init and script

