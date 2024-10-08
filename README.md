# PSM_Explants
This code simulates the growing PSM explants.

Agent-Based Model for Viscous Tissue Growth in a Porous Matrix:

This code, written, developed, and maintained by  Dr. Anupam Gupta, simulates the growth of embryos in a variety of confinement conditions using an agent-based model. It builds upon the work presented in https://doi.org/10.1242/dev.199423, by incorporating the possibility to add different confinement.

Features:

Modular Fortran Code: Written in Fortran for compatibility with various Fortran compilers.

User-Friendly Input: The para.in file allows you to specify parameters for simulating growing/non-growing embryos, using a switch.

Key Parameters:

tau (line 20): Slowest time scale, controls the decay of incoming motile cells.

FGF_grad (line 46): Control, whether there should be a gradient in the motility along A-P direction or not.

TB_Wall (line 53): Yes: Ensures that there are rigid walls around the explant. No: Soft boundary of beads around the explant/embryo.

Sample Usage:

The provided sample parameters in para.in are configured to simulate a soft boundary that confines an explant with decay in motility along A-P direction.
