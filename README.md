# Flying Touchdown 

Flying touchdown code developed by Siddhesh Sakhalkar (siddhesh_sakhalkar@berkeley.edu). see /docs/Flying Touchdown Code Manual.pdf for steps to use this code.

## Overview

The Flying Touchdown Code computes the slider temperature profile, disk temperature profile and the slider fly height for a flying TFC slider over a rotating AlMg PMR disk. Simulations can be performed over a range of TFC powers starting from 0 mW (TFC heater is not powered) all the way down to contact and beyond. To accurately predict the fly height and heat transfer in the head-disk interface (HDI) at near-contact, we incorporate the effects of disk temperature rise, asperity-based adhesion & contact force models, air & phonon conduction heat transfer and friction heating in our model.

The flying touchdown code is based on the CML TFC Code (https://github.com/cml-berkeley/tfc) . It uses a modified version of the CML Air static simulation program (quick code) to solve the to solve for the air bearing pressure distribution and the steady-state fly height. CML Air determines an equilibrium flying state (fly height, pitch, roll) of a slider for a given suspension load using a Quasi-Newton method. The air bearing pressure is computed by solving the generalized Reynolds equation with the Fukui–Kaneko slip correction. In the modified CML Air quick code, the equilibrium flying state of the slider is then determined by balancing the forces and torques on the slider due to adhesion and contact forces determined using the Sub-Boundary Lubrication model, the suspension load and the air bearing pressure.

The slider temperature profile and thermal protrusion is determined using a thermo-mechanical finite element (ANSYS) model, similar to the CML TFC code. The heat transfer coefficient due to pressurized air cooling is obtained using the temperature jump theory and the modified mean free path of air due to boundary scattering. The heat transfer coefficient due to wave-based phonon conduction is determined as a function of the spacing, the temperatures and material properties of the head and the disk. The net heat transfer coefficient in the HDI due to air & phonon conduction is integrated into the slider ANSYS model to simulate the head temperature profile and head thermal protrusion due to TFC heating.

The temperature profile of the rotating disk due to heat transfer from the head is determined using the analytical solution of the classic “stationary heat source acting on the surface of a moving semi-infinite medium” problem [9]. The net heat generation rate per unit area due to friction is determined as: 𝑞𝑓𝑟𝑖𝑐 =𝜇𝑝𝑈. Here 𝜇 is the coefficient of friction, 𝑝 is the net normal pressure due to the contact and adhesion forces and 𝑈 is the linear disk speed. For further details, please read Ref. [1].

The Flying Touchdown code runs on a windows-platform computer with MATLAB and certain ANSYS product (modules for electric, thermal and structural analysis) installed.

Ref 1: Sakhalkar, S., Cheng, Q., Ghafari, A. and Bogy, D., 2020. Investigation of heat transfer across a nanoscale air gap between a flying head and a rotating disk. Journal of Applied Physics, 128(8). (Available at [CML website](https://cml.berkeley.edu/cml-blue-reports/) as a blue report)

