# README: Flashing Atomization UDFs for ANSYS Fluent

## Overview
This repository contains User Defined Functions (UDFs) developed for ANSYS Fluent to simulate the **flashing phenomenon** of cryogenic propellants (liquid nitrogen and oxygen) under vacuum conditions. The model extends the work of Adachi et al. [1] and Zuo et al. [2] by incorporating radiative heat transfer, which is critical for spacecraft propulsion systems. The code enables simulations of superheated sprays, including droplet vaporization, heat/mass exchange, and spray pattern evolution.

## Key Features
- **Flashing Vaporization Model**: Accounts for internal superheat, conductive/convective heat transfer, and radiative heat exchange.
- **Euler-Lagrange Approach**: Tracks droplets in a Lagrangian frame while solving the continuous phase in Eulerian frame.
- **Validation**: Matches experimental spray patterns and temperature profiles from cryogenic propellant experiments.
- **Customizable**: Adjustable parameters for droplet size distribution, injection conditions, and turbulence modeling.


## References
1. Meng Luo and Usman Rana Modeling Investigation of Liquid Oxygen Flashing Spray with CFD
2. Thesis: *Numerical Investigation of Cryogenic Propellants’ Flashing Phenomenon* (Usman Rana, 2017).
1. Adachi et al. (1997), *SAE Technical Paper*.  
2. Zuo et al. (2000), *International Journal of Engine Research*.  

## License
This code is part of a thesis project. All rights reserved. For permissions, contact the author or supervising institution.  

---

**Author**: Usman Rana  
**Contact**: ur19@ic.ac.uk  
**Thesis Supervisor**: Meng Luo, M.Sc  
**Institution**: Technische Universität München  
