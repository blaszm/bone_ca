# bone_ca
Program to perform multiphysics, multiscale bone simulations for research purposes

## Description
Current version updated on: 27th February 2024

Osteoporosis is the leading bone disease worldwide and especially prevalent in the elderly population. A possible future method of early detection could be the usage of sonography. Our research focuses on multiphysiscs, multiscale simulations of human bone, which could aid the development of new diagnosis tools.

In our work, cancellous bone is modeled as a composite material consisting of the phases cortical bone and bone marrow. A multiscale approach (FE2) is used in conjuction with representative volume elements (RVEs) on the microscale to model different stages of osteoporosis and to compare the results. On the macroscale, a boundary box filled with air is included to account for a surrounding medium and to enable the decay of the resulting electro-magnetic fields.   

## Usage
For microscale simulations, run the main document or bone_microscale_main.jl. The example compares RVEs, computing the effective Young's modulus and the average flux quantities for healthy and ill bone.
The macroscale simulation can be run by executing bone_macroscale_main.jl (very long calculation, which is parallelized, **requires a computer cluster**).

## Journal publications

1. [Multiscale modeling of cancellous bone considering full coupling of mechanical, electric and magnetic effects](https://doi.org/10.1007/s10237-021-01525-6)
2. [Inverse modeling of cancellous bone using artificial neural networks](https://doi.org/10.1002/zamm.202100541)
3. [On the effects of a surrounding medium and phase split in coupled bone simulations](https://doi.org/10.1002/zamm.202200595)

For further reading and information, visit my [Orcid](https://orcid.org/0000-0002-1126-3303).
