# Shallow Water Simulation Unity Project

## Project Overview

This is a Unity project for shallow water simulation and rendering, using GPU compute shaders to accelerate numerical calculations.  
Key features include:

- **GPU-accelerated simulation**: Compute water height and velocity in parallel on the GPU.  
- **Flexible parameters**: Adjustable water depth, gravity, grid size, grid resolution, and other physical parameters.
- **Wet/dry supported**: A rough but stable dry–wet treatment is implemented. 
- **Efficient grid handling**: The simulation is stored in 1D arrays while mapping clearly to 2D grids, and enabling straightforward half-grid flux calculations.

---

## Demo Video

Bilibili: https://www.bilibili.com/video/BV1ZoFVzDE2Q/

YouTuBe
- Waves on the beach: https://youtu.be/VkTo75j1LR8
- Waves with Hills: https://youtu.be/4MepGrY8EnU

---

## Project Structure

```
  ShallowWater/
  │
  ├─ Assets/ # Main Unity assets
  │   ├─ Scripts/ # C# scripts
  │   │   └──ShallowWaterCPUControl.cs # Main simulation control, buffer management
  │   └─ Shaders/ # Compute and rendering shaders
  │       ├──ShallowWaterSim.compute # Main shallow water solver (HLLC Riemann sovler, Forward Euler method)
  │       ├──WaterHeight.shader # Water surface render shader
  │       └──SeaBed.shader # Seabed render shader
  └─ README.md # Project description
```

---

## Author

Developed by **Shunquan Huang**  
Ph.D. Candidate in Astrophysics, UNLV  

---

## License

This project is released under the **MIT License**.  
Feel free to use, modify, and distribute with attribution.
