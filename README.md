# Configuration Design of an EDF-Driven Personal Flight Suit: Optimizing for Power and Aerodynamic Efficiency
_Alecia Barbieri — Princeton University Senior Thesis, 2025_

## Project Overview

This repository contains the final design of my flight suit and supporting materials for my senior thesis in Electrical and Computer Engineering at Princeton University. The project explores the design of a compact, efficient, and lightweight ducted electric propulsion (EDF) glide system to enable low-speed personal flight.

The thesis integrates principles of aerodynamics, propulsion sizing, and power distribution to produce a system that is aerodynamically optimized for the human form and capable of generating sufficient thrust for a 120 lb person during hover and cruise flight.

## Objectives
- Iteraively optimize aerodynamic parameters for maximized lift-to-drag ratio while minimizing weight and wingspan for human-interface ergonomic design.
- Utilize XRotor data sweeps to maximize propeller performance and match to motor with characteristic yielding maximized efficiency across hover, climb, and cruise flight profiles.
- Model power consumption and mission energy needs for different flight phases and sizing power source to achieve those requirements.
- Estimate structural mass and select optimal materials for fixed wing and EDF components.
- Validate performance through theoretical modeling and numerical tools.

## Repository Contents
- **PersonalFlightSuitDesigns Folder:** Contains the different iterations of propeller geometries, including distributed lift coefficient and twist (.txt files can be loaded into XROTOR with the load command) and finalized OpenVSP geometry, including the human, pod-man, finalized wing geometry (BRC42), and the finalized EDF design (Prop20)
- **PolarData Folder:** Contains VSPAero polar arrays for each wing iteration
- **XROTORSweepData Folder:** Contains RPM, velocity, and advanced ratio sweeps from propeller design, extracted from XRotor
- **ThesisOptimizationVariables.xlsx:** An Excel notebook that contains all of the saved iterations of the wing geometry, propeller geometry, and considered motors
- **DuctedPotentialFormulationSolutions.txt:** Contains miscallaneous data on each propeller design
- **MATLAB scripts:** Used to analyze data and generate plots to visualize and validate my calculations

## Getting Started

To replicate or extend the analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/aleciabelle6/AleciaBarbieriPrincetonSeniorThesis.git
   cd AleciaBarbieriPrincetonSeniorThesis
   ```

2. Open MATLAB and add this new folder to your path.

## Requirements

- MATLAB R2023b or later
- OpenVSP (and in-application VSPAero)
- XRotor

## MATLAB Code
Scripts are executed in the following order, with any subsequent script running all of the scripts before it:
1. FlightEnvelope.m - Initial Mass Estimates and Wing Geometry Optimization (Helper Script: Glide.m)
2. Mission.m - Defines the mission, specifically time intervals and distances
3. Propulsion.m - Determines Thrust Required for each Segment: Cruise, Optimal Glide, Hover
4. PropDesign.m - Sizing and Optimizing Propeller for Mission Requirements (Helper Script: MotorSelection.m)
5. PowerRequired.m - Calculates Power Requirements for our sizing our Power Source (Helper Script: TransitionPeriod.m)

PropDesignCompareProps.m was developed after the submission of my report to automate both the propeller and motor efficiency comparison (as opposed to just pooling motors like in PropDesign.m) and, when ran, also runs files 1-5 prior to executing in full.

## Wing Optimization Iteration Process
To run the wing optimization algorithm:
- Create a wing design in OpenVSP
- Extract the polar data to retrive the aerodynamic coefficients (stored in the PolarData folder)
- Add another column to the main ThesisWingOptimizationVariables.xlsx and manually input the geometry parameters (light and dark orange, except for density and mass)
- Extract the density and input into OpenVSP and run the Mass Calculator to extract the mass (update spreadsheet)
- Run FlightEnvelope.m with the variable _iteration_ being the prop number of your current iteration to autopopulate the rest of the variables into the spreadsheet

## Propeller Optimization Iteration Process
To run the propeller optimization algorithm:
- Design a propeller informed by the thrust requirements from Propulsion.m
- Extract the RPM, thrust, power, and efficiency sweeps for the 3 different mission profiles (stored in the XRotorSweepData folder)
- Add data from XRotor Main menu for that design to the bottom of DuctedPotentialFormulationSolutions.txt
- Enable Propeller Geometry customization on FlightEnvelope.m in determining initial weights and set _iteration_ to the chosen wing you want to run with this propeller.
- PropDesign.m with the variable _propIteration_ being the prop number of your current iteration to match that propeller to pooled motors to optimize for efficiency

## Example Figures

Here are some key outputs from the analysis:

- W/P vs. W/S Flight Envelope
- C_D vs. C_L Polar Curve for Calculating Optimal Glide Angle
- Comparison of Aerodynamic Parameters employing Panel vs. VLM method on VSPAero
- Motor Efficiency vs RPM
- Power Required vs. RPM
- Power Required vs. Number of EDFs on Flight Suit

## License

This project is made available for academic, non-commercial use. Please contact me if you'd like to use any part of the content or code beyond that scope.

---
*Author: Alecia Barbieri ([LinkedIn](https://www.linkedin.com/in/aleciabarbieri))*  
*Project Advisor: Dr. Luigi Martinelli*  
*Department of Electrical and Computer Engineering, Princeton University*
