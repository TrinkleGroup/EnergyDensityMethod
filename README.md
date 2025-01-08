<h1 style="text-align: center;"> Energy Density Method (EDM) Implemented in Vienna Ab initio Simulation Package (VASP)</h1>

The Energy Density Method (EDM) is a theory based on density-functional theory (DFT) that decomposes the Kohn-Sham total energy of a solid-state system into well-defined atomic energies. This energy decomposition involves two steps described as below.

* **Step 1:** Calculate an energy density function $e(\mathbf{r})$ based on DFT.
  - The energy density function restores the Kohn-Sham total energy, up to a constant difference, if integrated over the entire real space, $E_\text{Kohn-Sham} = \int d\mathbf{r} e(\mathbf{r}) + C$.
  - Also, the potentials and charge densities useful for step 2 are calculated and exported.
* **Step 2:** Integrate the energy density components over volumes partitioned to individual atoms to obtain atomic energies.
  - Certain conditions in the volume partitionings should be met to ensure the **uniqueness** (gauge-invariance) of the atomic energies.
  - These conditions are met by introducing **Bader volumes** and **charge-neutral volumes**, based on Bader analysis of the charge density and the total potential given by step 1.
  - The atomic energies $E_{\text{EDM}, \mu}$, if summed over all atoms, restores the Kohn-Sham total energy up to a constant difference, $E_{\text{Kohn-Sham}} = \sum_\mu E_{\text{EDM}, \mu} + C$.

This code performs the first step, completing standard DFT calculations set up by VASP while doing some extra work to obtain the energy densities, charge densities and potentials. For the second step, a Bader analysis/integration code is needed. You can find a usable one in the [Bader Integration Repository](https://github.com/TrinkleGroup/BaderIntegration) developed by the [Trinkle Research Group](https://github.com/TrinkleGroup).

Please see the [Reference](#references) part for detailed explanation and mathematical derivations of this method.

Below is a guideline for users who wish to perform EDM calculations.

# Compilation
EDM is implemented in the framework of [Vienna Ab initio Simulation Package (VASP)](https://www.vasp.at/). VASP is a commercial software and a valid license is required for authorized uses. The implementation of EDM is presented as a patch to the original VASP source code. The current version of EDM code is based on VASP version 5.4.4. The compilation instructions are as below.

1. Prepare the VASP source code and check its installation requirements in your machine. For the source code, make sure that all the bug-fixing patches have been applied. VASP 5.4.4 needs to be patched with file ```patch.5.4.4.16052018```. See the [VASP manual for installation](https://www.vasp.at/wiki/index.php/Installing_VASP.5.X.X) for details.

2. Place the EDM patch file, for example, file ```EDM_VASP.5.4.4_v3.2_stable.patch``` in the root directory of VASP's source code. Run

       patch -p1 < EDM_VASP.5.4.4_v3.2_stable.patch

   to get the EDM source code.

3. Prepare file ```makefile.include```. Take the file ```makefile.include``` that compiles pristine VASP in your machine, and add the following precompiler flags to the ```CPP_OPTIONS```: ```-DSUPPORT_EDM```, ```-DSUPPORT_EDM_paw```, ```-DSUPPORT_EDM_GGA```, so that it looks something like

       CPP_OPTIONS= ... \
                    -DSUPPORT_EDM \
                    -DSUPPORT_EDM_paw \
                    -DSUPPORT_EDM_GGA

   These flags control the precompilations of the EDM code in general, the EDM code containing projector augmented-wave (PAW) method, and the EDM code containing the generalized gradient approximation (GGA). 

   Note:
      * All EDM-related operations are coded under statement ```#ifdef SUPPORT_EDM```, so if flag ```-DSUPPORT_EDM``` is missing, the compiler will not compile any EDM code and will do a normal VASP compilation.
      * All PAW- and GGA-related EDM code parts, controlled by ```-DSUPPORT_EDM_paw``` and ```-DSUPPORT_EDM_GGA```, respectively, are nested inside ```-DSUPPORT_EDM```. They exist primarily for development and for historic reasons and are currently kept in the code.

4. Compile the EDM code. In the root directory of the source code, run

       make std

   to make the binary file ```vasp_std``` under the ```root/bin``` directory, which will be the binary file for EDM calculations. Note that EDM is not implemented in gamma-only, non-collinear or GPU versions of VASP.

# Input
In general, the input files are the same as what are needed by VASP calculations: INCAR, KPOINTS, POSCAR and POTCAR, so please follow the guidelines in [the VASP manual](https://www.vasp.at/wiki/index.php/The_VASP_Manual) to prepare the inputs. However, additional requirements are needed for the INCAR file to perform accurate and efficient EDM calculations, which are listed as below.

## Required (should be applied at all times or EDM will encounter problems):
1. Manually set the FFT-grid and the "fine" FFT-grid to equal numbers of grid points in each dimention, i.e., ```NGX=NGXF```, ```NGY=NGYF```, ```NGZ=NGZF```. 
   * This is required for successful energy integration for atomic energies.
   * EDM requires denser real-space grids than typical DFT calculations to achieve good accuracy of atomic energies, so please set these values to be higher than the default VASP recommendations. For example, with ```PREC=Accurate```, by default VASP sets ```NGX=2Gcut``` and ```NGXF=2NGX```, etc. A common practice is to set ```NGX=NGXF=4Gcut```, etc. Same for ```NGY(F)``` and ```NGZ(F)```.
2. ```LELF = .TRUE.```
   * Kinetic energy densities and their gauge-dependent terms are updated in subroutine ELF, which is not enabled unless ```LELF=.TRUE.``` is turned on.
   * VASP requires that ```NPAR=1``` is set with ```LELF=.TRUE.```. However, we've allowed ```NPAR>1``` for more efficient parallelized EDM calculations, and made sure that the EDM calculations are correctly parallelized; avoid using the ```ELFCAR``` file in this case.
3. ```NSW = 0``` 
   * One should not perform any relaxations with EDM. EDM will not raise warnings or errors if ```NSW>0``` but the energy density data are not reliable in this case. Also, as a side note, one should avoid using the energy data calculated in ionic relaxations by VASP for purposes that require accurate energy of a geometry; in such cases, static calculations should be performed.
4. ```IBRION = -1``` 
   * Required with ```NSW=0```, so that no ionic updates are performed.
5. ```ICHARG != 4```
   * Setting ICHARG=4 will bypass the code in which EDM calculations are done. 
   * EDM is not implemented for OEP.
    
   ### And the default settings for the following tags should be kept.
6. ```KPAR = 1``` 
   * EDM is not implemented nor well-tested for k-mesh parallelism with ```KPAR>1```.
7. ```IVDW = 0```
   * EDM is not implemented nor well-tested with vdW correction.
8. Possibly others. If a specific feature or correction needs to be manually turned on in ```INCAR``` and not mentioned in this EDM manual, it is likely that EDM is not implemented for that feature/correction, and may or may not capture its energy effects (if any) in the energy densities. Please read the source code, do your own test, or contact us if you are not sure.

## Recommended (should be applied under certain conditions for better performance):

1. ```ADDGRID = .TRUE.```
   * We recommend using it when possible for better accuracy, as EDM usually needs denser real-space grids than what VASP recommends by default.
   * However, be careful of aliasing errors. User should check the [VASP Wiki on ADDGRID](https://www.vasp.at/wiki/index.php/ADDGRID) for suggestions and warnings of potential problems.

2. ```LREAL = .FALSE.```
   * We recommend using projection in reciprocal space in EDM calculations for better quality in energy data whenever possible.
   * Real-space projections always yield an error in the DFT energy, which is small but not necessarily negligible, and is usually a constant energy shift for each atom. If ```LREAL``` is not set to ```.FALSE.```, the user **should** make sure that a proper reference EDM calculation is done with the same ```LREAL``` setting, among other settings. See [VASP Wiki on LREAL](https://www.vasp.at/wiki/index.php/LREAL) for more information.
   * It is possible to apply ```LREAL=.FALSE.``` even to large systems in VASP or EDM calculations. For example, one can perform the geometry relaxation and electronic optimization using ```LREAL=Auto```, and then do a final electronic optimization and EDM calculation using ```LREAL=.FALSE.```. See the [Assorted Wisdom](https://rosengroup.slite.page/p/PjutyhI3iHVl7c/Assorted-Wisdom) and [Input Files](https://rosengroup.slite.page/p/0CKuG1eLlFb3Lc/Input-Files) from [A Practical Guide for VASP, Rosen Review](https://rosenreview.cbe.princeton.edu/dft/vasp) for more information and comments on ```LREAL```.
3. ```LAECHG = .TRUE.```
   * If the PAW method is used, we recommend generating all-electron charge densities, which could be useful in further analyses of gauge-invariant volumes and EDM atomic energies.

## Outdated requirements and recommendations (required/recommended for older versions of EDM only)
1. ```NPAR = 1```
   * Required for EDM implemented in VASP 4.6, and for EDM implemented in VASP 5.4.4 prior to versions series v2.
2. ```ALGO = Normal``` 
   * Required for EDM implemented in VASP 4.6, and for EDM implemented in VASP 5.4.4 prior to versions series v3.
3. ```LMAXMIX = 4```
   * Recommended for EDM implemented in VASP 4.6, and for EDM implemented in VASP 5.4.4 prior to versions series v3, for convergence stability with *d* or *f* orbitals. EDM has been restructured to be a pure post-processing operation starting from version series v3 and is no longer involved in the SCF loops.

# Output
The output files contain the standard VASP outputs as well as the EDM-specific files listed below. The EDM-specific files have a similar format as the nonmagnetic [CHGCAR file](https://www.vasp.at/wiki/index.php/CHGCAR), containing
* Structural information in [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) format;
* FFT Grid dimensions, ```NGX*NGY*NGZ``` or ```NGXF*NGYF*NGZF```;
* FFT Grid volume times the EDM-generated quantity (energy density, charge density or potential) on the grid.

For spin-polarized calculations, the spin-polarized quantities are appended with "\_1" and "\_2" in the filenames, denoting spin-up and spin-down quantities. The total quantities are appended with "\_tot" in filenames, which are the direct sum of the corresponding "\_1" and "\_2" quantities.

## Energy density files

1. ```EDM_Ta```: kinetic energy density (the asymmetric form), $t^{(a)}(\mathbf{r})$.

2. ```EDM_Tc```: gauge-dependent term of the kinetic energy density, which is the difference between the asymmetric and the symmetric form, $t^{(a)}(\mathbf{r}) - t^{(s)}(\mathbf{r}) = -\frac14 \nabla^2 \tilde{\rho}(\mathbf{r})$.

3. ```EDM_VnCAR```: total classical Coulomb energy density (the asymmetric form), $e_\text{CC}^{(a)}(\mathbf{r})$.

4. ```EDM_exc```: exchange-correlation energy density, $e_\text{XC}(\mathbf{r})$.

5. ```EDM_ONSITE```: non-local on-site energies assigned to each atom $\mu$, $E_\mu^\text{on-site}$. For spin-polarized calculations, the file contains two columns, corresponding to the spin-up and spin-down energy components, and the total atomic on-site energies are the sum of these two columns.
    
## Charge density files

1. ```EDM_CHDEN```: soft pseudo-charge density $\tilde{\rho}(\mathbf{r})$, used to define Bader volumes.

2. ```EDM_NTCAR```: the charge density part in the equation for classical Coulomb energy density, $\rho^e(\mathbf{r}) - \rho^\text{model}(\mathbf{r})$.

## Potential files

1. ```EDM_VCCCAR```: the total classical Coulomb potential $V^\text{tot}(\mathbf{r})$ times -1, $-V^\text{tot}(\mathbf{r})=-\left[V_\text{H}(\mathbf{r})+V^\text{loc}(\mathbf{r})\right]$, used to define charge neutral volumes.

2. ```EDM_VTCAR```: the potential part in the equation for classical Coulomb energy density, $2V^\text{loc}(\mathbf{r}) + V_\text{H}(\mathbf{r}) + V^\text{model}(\mathbf{r})$.


## OUTCAR
The summed and integrated kinetic energy $T$, classical Coulomb energy $E_\text{CC}$, exchange-correlation energy $E_\text{XC}$ and non-local on-site energy $E_\mu^\text{on-site}$ are written near the end of the ```OUTCAR``` file. The summations are done over all atoms $\mu$ and the integrations are done over the entire real space, so these quantities reflect the _total_ EDM energy components, which sum up to be the EDM total energy. For spin-polarized calculations, spin-polarized energy components are also presented. The EDM total energy differs from the Kohn-Sham total energy by a constant from the pseudopotential, which equals $\sum_\mu \int d\mathbf{r} V_\mu^\text{loc}(\mathbf{r}) \rho_\mu^\text{loc}(\mathbf{r})$.

## Verbose-mode files
Verbose mode (the previous debug mode) can be activated by setting ```EDMVERBOSE=1``` in ```INCAR```, which will generate some extra files for EDM quantities.

1. ```EDM_VloCAR```: the local part of the pseudopotential, $V^\text{loc}(\mathbf{r})$.

2. ```EDM_VhCAR```: Hartree potential, $V_\text{H}(\mathbf{r})$.

3. ```EDM_VmodCAR```: model potential, $V^\text{model}(\mathbf{r})$.

4. ```EDM_nmodCAR```: model charge density, $\rho^\text{model}(\mathbf{r})$.

5. ```EDM_neCAR```: valence electron density $\rho^e(\mathbf{r})$.


# Typical workflow for EDM calculations and analyses
EDM should be considered as a post-processing of the outputs and parameters from the DFT calculation. In other words, standard DFT calculations should be performed to relax the geometry and converge the wavefunction and charge density before EDM takes effect. In the most recent version, EDM calculations and memory allocations are performed only after the electronic optimization is complete. A typical workflow of EDM is listed as follows.

1. Geometry optimization, using standard DFT.

2. Self-consistency calculation (```NSW=0```, ```IBRION=-1``` and ```NELM>0```) for converged wavefunction and charge density, using standard DFT. 
   * Note that most other settings should be the same with the upcoming EDM calculation for consistency.

3. Take the converged wavefunction and charge density as input, restart (```ISTART=1```) an EDM calculation with electronic optimization turned off (```NELM=0```, which means that the wavefunction and charge density are kept frozen). Check the [Input](#input) part for guidance on the ```INCAR``` settings.
   * This step can actually be combined with the previous step by running EDM with ```NELM>0```. However, for very large systems, we recommend considering to perform these two steps separately. The reason is that there is a surge in memory use when the EDM calculation begins, which, if the allocated job memory is already tight, could break the job down before the results are properly saved, after a costly run to reach SCF convergence.
   * The converged charge density, if absent, can alternatively be constructed from the converged wavefunction.

4. Numerically integrate the energy densities over gauge-invariant volumes for atomic energy components. Note that each energy density should be matched with the correct gauge-invariant volumes: Bader volumes (found through the pseudo charge density in file ```EDM_CHDEN```) for the kinetic and exchange-correlation energy densites, and charge-neutral volumes (found through the total classical Coulomb potential in file ```EDM_VCCCAR```) for the classical Coulomb energy density. If using the code for [Bader integration with weight method](https://github.com/TrinkleGroup/BaderIntegration), the commands are

       weight_int -s EDM_CHDEN EDM_Ta > ta.txt
       weight_int -s EDM_VCCCAR EDM_VnCAR > ecc.txt
       weight_int -s EDM_CHDEN EDM_exc > exc.txt

   which generates the kinetic, classical Coulomb and exchange-correlation atomic energies, respectively. Similarly, the integration error can be assessed by integrating the gauge-dependent term via command

       weight_int -s EDM_CHDEN EDM_Tc > tc.txt
    
    For spin-polarized calculations, the atomic energy components are found by running

       weight_int -s EDM_CHDEN_tot EDM_Ta_tot > ta_tot.txt
       weight_int -s EDM_VCCCAR EDM_VnCAR > ecc.txt
       weight_int -s EDM_CHDEN_tot EDM_exc_tot > exc_tot.txt

    (For checking purposes only in most cases) Spin-polarized atomic energy components are obtained with commands like

       weight_int -s EDM_CHDEN_1 EDM_Ta_1 > ta_1.txt

5. Perform a reference EDM calculation for reference energy, with settings consistent with steps 1-4.

   The absolute values of raw EDM atomic energy data lack clear physical meanings, and for practical use, a reference energy from EDM calculation on a properly chosen reference system is *always* necessary.

6. Subtract the reference energy components obtained from step 5 from the atomic energy components obtained from step 4, and sum up the energy components to get the EDM atomic energies.

# Version history
## Implementation in VASP 5.4.4
### V3 series
* v3.4 (beta version)
   - Bug fix for the meta-GGA part in v3.3 that potentially breaks down spin-polarized EDM calculation with normal GGA.
* v3.3 (obsolete)
   - Added functionality for nonspin-polarized EDM calculation with meta-GGA.
* v3.2 (stable version)
   - Added improvement over numerical stabilities and file IOs, mostly for classical Coulomb energy calculations; problem fixes.
* v3.1 (obsolete)
   - Code structure reorganized and optimized for improvement in efficiency in CPU time and memory, primarily for EDM calculations on very large (>500 atoms) systems, such as dislocations.
      - EDM operations were entirely moved out of the electronic optimization loop, making it possible to perform EDM calculations with frozen wavefunctions and charge densities, and removing the constraint on the ```ALGO``` tag.
      - EDM operations and file IOs were minimized, so that each EDM quantity is updated only once, and each EDM file is read (if necessary) and written only once.
      - Memory was optimized to reduce the maximum required memory.
      - Full parallelization was enabled over electronic bands, at the same level as the original VASP code.
   - On the largest calculation (917 atoms or 7336 valence electrons, spin-polarized) that we've tested, this update reduced the runtime from weeks (with V2 version series of EDM) to around 12 hours.

### V2 series (obsolete)
* Initial implementations and full tests of EDM with (collinear) spin polarization (```ISPIN=2```);
* Code structures similar to v1 series;
* Constraint of ```NPAR=1``` removed for better parallelization over electronic bands.

### V1 series (obsolete)
* Initial implementations of EDM in VASP 5.4.4, with similar code structures, functionality and input requirements compared to the implementation in VASP 4.6;
* Nonspin-polarized EDM with ```ISPIN=1```;
* Parallelization over electronic bands with ```NPAR=1```;
* EDM data are updated in the electronic optimization loop with the blocked-davidson method (```ALGO=Normal```).

## Implementation in VASP 4.6
### V0.0 (archived)
* Initial implementation of EDM by Min Yu in VASP 4.6.36 and VASP 4.6.38.

# References
* M. Yu, D. R. Trinkle, and R. M. Martin, "Energy density in density functional theory: Application to crystalline defects and surfaces." *Phys. Rev. B* **83**, 115113 (2011). [doi](http://dx.doi.org/10.1103/PhysRevB.83.115113)
* More references to come.

# Contributors
Yang Dan (version 5.4.4), Min Yu (version 4.6.36/4.6.38) and Dallas R. Trinkle