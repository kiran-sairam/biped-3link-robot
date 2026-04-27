# 3-Link Planar Bipedal Robot Simulation & Control (HZD)

This repository contains the modeling, optimization, and simulation MATLAB files for the three link bipedal robot, for the mini project requirements for the skpring 2026 Legged Robotics course.

The locomotion control strategy uses **Hybrid Zero Dynamics (HZD)** and **Feedback Linearization**, using Bézier polynomials to enforce virtual constraints for stable, periodic walking gaits.

## Overall Mapping

The codebase is organized to show the four projects:

### Mini-Project 1 & 2: Dynamic Modeling & State Space Representation
**Aim:** To derive the forward kinematics, Lagrangian dynamics, and state space representation for the single support phase and impact mapping.
* **`genrate_functions.m`**: The symbolic file that derives the kinetic/potential energies, the $D$, $C$, $G$, and $B$ matrices, and the extended impact mapping Jacobians.
* **`autogen/`**: Contains the optimized, auto generated MATLAB functions (e.g., `func_compute_D_C_G_B.m`, `func_compute_De_E_dY_dq.m`) exported directly from the symbolic derivations.
* **`func_full_dynamics.m`**: Puts together the state space formulation $\dot{x} = f(x) + g(x)u$.

### Mini-Project 3: Zero Dynamics & Gait Optimization
**Aim:** Partitioning the dynamics to isolate the zero dynamics equations and formulate an optimization problem to find the most energy efficient walking gait.
* **`func_zero_dynamics.m`**: Evaluates the 2D zero dynamics manifold using dynamic partitioning.
* **`Optimize.m`**: The nonlinear programming script (`fmincon`). Optimizes the pre impact states and Bézier coefficients to minimize the Mechanical Cost of Transport (MCOT) while staying within the ground constraints.
* **`sim_zero_dynamics.m` & `sim_and_plot_ZD.m`**: Simulates the reduced order zero dynamics to verify limit cycle existence before using full body control.

### Mini-Project 4: Nonlinear Control & Stability Analysis
**Aim:** Design a feedback linearizing controller to move the system to the zero dynamics manifold and verify stability using Poincaré phase portrait method.
* **`func_feedback.m`**: Implements output feedback linearization coupled with a PD controller.
* **`func_impact_map.m`**: Calculates the instantaneous rigid body impact dynamics and relabels the state coordinates for next steps.
* **`sim_and_plot_full_dynamics.m`**: The primary execution script. Integrates the continuous swing phase and discrete impacts over multiple steps to validate the controller.
* **`plot_trajectories.m` & `animate_results.m`**: Generates the phase portraits ($q_1$ vs. $\dot{q}_1$) and renders the physical animation of the biped.

---

## Quick Start Execution

To run the full simulation and generate the stability phase portraits:
1. Run `set_path.m` to add the utility and autogen directories to your path.
2. Rename filepaths to your local PC's MATLAB directory.
