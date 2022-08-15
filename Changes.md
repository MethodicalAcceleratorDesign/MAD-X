MAD-X master

*   [PR 1125](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1125) Make TWISS treat [XY]ROTATION exactly, including the linear and second order map (J. S. Berg)
*   [PR 1123](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1123) Fix errors in the time variable with EXACT flag to TWISS (J. S. Berg)
*   [PR 1081](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1108) Fix potential buffer overrun when node_name calls stoupper (J. S. Berg)
*   [PR 1107](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1107) Stabizes few tests due to compiler dependent numerical noise (R. De Maria)
*   [PR 1095](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1095) Implement more robust, optional, PTC DA map output (L. Deniau)
*   [PR 1088](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1088) Additional explanation bv flag (J. Dilly)
*   [PR 1113](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1113) Introduce LAST option in INSTALL and MOVE and change default behaviour [Breaking change!]
*   [PR 1131](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1131) Add spin table (P. Skowonronski)
*   [PR 1093](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1093) Refine aperture refinement(T. Persson)
*   [PR 1114](https://github.com/MethodicalAcceleratorDesign/MAD-X/pull/1114) Exact option in translation(J. S. Berg)

MAD-X release 5.08.01 (2022.02.25)

*    Fix radiation of mulitpoles in TRACK. #1079. (Riccardo)
*    Fix in exact drift transfer map in TWISS #1077. (Tobias)
*    Allows disabling the scaling of TWISS in PTC by default but possible to activate it. #1073. (Tobias)
*    Fixes a bug in the dqmin calculation. #1075. (Tobias) 

MAD-X release 5.08.00 (2022.01.17)

*    Wire is now fully implemented in TWISS and TRACK. #1018. (Tobias, testing: Phillipe and Guido)
*    Bug in the output of fc.3 where sometimes the highest orders multipoles were not included. #1022. (Riccardo)
*    Restore the internal names in the SURVEY table. #1024. (Tobias)
*    When an element was replaced the keyword was not updated. This is fixed now. #1024. (Tobias)
*    Bug in TWISS with centre option combined with a tilted quadrupole. #1026. (Tobias)
*    Bug in EMIT for octupole where k2 was retrived instead of k3. #1028. (Ghislain)
*    Exact drift option implemented in TWISS. #1030. (Tobias)
*    New release of PTC/FPP. #1031. (Piotr, Tobias)
*    Npart for beam-beam element can now also be set directly in the element. #1033. (Tobias)
*    New implementation of the tapering. #1034. (Ghislain, Tobias)
*    Array potentially out of range. #1041. (J.S. Berg)
*    Fix issue in tracking when FFILE is used. #1048. (Tobias)
*    The unified atomic mass is introduced (AMASS). #1049. (Tobias)
*    The information from EMIT is now stored in a table. #1053. (Tobias)
*    The solenoid was wrong with a factor 2 in EMIT. #1061. (Andrea, Riccardo)
*    Estimate of the DeltaQmin is now included in TWISS. #1060. (Tobias, Eirik)
*    Option to output TPSA from PTC/FPP with more digits and cleaner format. #1052. (Laurent, Tobias) 

MAD-X release 5.07.00 (2021.05.03)

*    Enable to save and load sequences with double precision floating point numbers in a hexadecimal format. In this way the exact same values are obtained after a saving and then loading a sequence. #954. (Tobias)
*    Translations and rotations are transported to PTC. #957. (Tobias, Laurent, J.S. Berg)
*    Voltage from crab cavities were previously included in total voltage exported to SixTrack. #965. (Tobias)
*    Possible to load a sequence and plot without any active sequence. #967. (Tobias)
*    Use Make 4.3 to build MAD-X. #969. (Laurent)
*    Too long sequence names got corrupted before now it is prevented. #971. (Tobias)
*    Lost particles in PTC caused wrong indexing resulting in wrong coordinates for lost particles (when more than 1). #972. (Tobias, Piotr)
*    Crabcavity with 0 length were ignored by MAD-X. #986. (Tobias)
*    Inconsitent use of beta and beta0 in TRACK. #988. (Andrea)
*    The tilt is now consistent in PTC with MAD-X and the documentation. #991. (Tobias)
*    Change to the geometrical definition for the Y-rotation in Survey (rest was already OK). #993. (Tobias, Laurent)
*    Introduction of permanent misalignments in MAD-X and a new survey option. #973. (Tobias)
*    Improved step size to find closed orbit in TRACK. #990. (Tobias)
*    Long names and more digits are now default in the SixTrack conversion. #999. (Tobias)
*    A bug with makethin with solenoid and expression was fixed. #1002. (Helmut)
*    Tapering of the phase removed, added to matching and a warning for unstable closed orbit in TWISS is added. #1005. (Tobias)
*    Correctors with length=0 but lrad will now radiate in TRACK. #1006. (Tobias)
*    Allowing the 3 parameter to be 0 in the definition of the RACETRACK. #1007. (Tobias)
*    Possibility to define a different beam-size-beam for horizontal and vertical in the APERTURE module. #1010. (Tobias) 

MAD-X release 5.06.01 (2020.09.01)

*    The BV flag was ignored for dipole fringe fields. This caused problems for LHC Beam2. Bug introduced in 5.06.00 release. #950. (Tobias) 

MAD-X release 5.06.00 (2020.08.14)

Added features/updates/improvements

*    Possible to write the eigenvectors of the one turn map to a file. #804. (Tobias)
*    Arbitrary aperture definition in TRACK. #806. (Tobias)
*    Signficant speed up of the thin TRACK. #807. (Tobias)
*    OPENMP compilation restored. #815. (Tobias, Laurent)
*    Signficant speed up of the thick TRACK. #817. (Tobias)
*    Implement flag to save the calculated corrector strength in sequence. #820. (Tobias)
*    Implement a new add expression command. #822. (Tobias)
*    Additional format to define arbitrary aperture geometry. #823. (Tobias)
*    When k0 is not equal to the angle is now exported to SixTrack. #826. (Tobias)
*    Thin combine function magnet implemented in TRACK and TWISS. #827. (Tobias, Malte)
*    Bug fix for the syntax of the option using ADD errors. #798. (Laurent)
*    Fix in PTC-MADX interface when a skew component was present togher with a normal component. #798. (Laurent)
*    3rd clock in PTC which enables RF-modulation. #833. (Piotr)
*    Arbitrary aperture geometry in PTC. #834. (Piotr)
*    Tilted aperture exported to SixTrack. #836. (Tobias)
*    Enables the plotting of longer sequences. #838. (Tobias)
*    If strength is 0 of a dipole it uses the exact drift in TRACK. #846. (Tobias)
*    PTC_TWISS with CENTRE option gave one element offset for the S postion. #859. (Piotr)
*    Get the average value when tracking many particles. #864. (Tobias)
*    PTC from 2020.03.25 #885. (Piotr)
*    When the CHROM option is used in TWISS the chromaticity is given as DQ/DPT (before it was DELTAP). #887. (Tobias)
*    Set Jacobian bisec default to 3. #893. (Laurent)
*    Fix of the NPART in BEAM that was overwritten when the beam was updated. #901. (Tobias)
*    RF-fringe field implemented. #902. (Tobias)
*    Tapering included in TWISS. #905. (Tobias, Leon)
*    Added WIRE attributes to the collimator and export to SixTrack. #914,#929. (Tobias, Helmut)
*    FDSTEP is now used in matching to define step size. #921. (Laurent)
*    Added radiation integrals I6 and I8. #924. (Tobias, Leon)
*    Fix warnings in plot module for non-twiss tables with haxis=s. #923. (Laurent)
*    Update constants to CODATA 2018. #925. (Tobias)
*    The strength of the dipole, if defined, is now defining the fringe fields for dipoles. #934. (Tobias) 

Bug Fixes

*    Fixed a bug in the move command when FROM was used for location. #815. (Tobias)
*    Segmentation fault when weights were used wrongly with MACRO. #835. (Tobias)
*    Fixes a bug in the MOVE command for complicated references. #854. (Tobias)
*    Radiation for a thin solenoid gave NaN values in TWISS. #858. (Tobias)
*    Fix of a problem with dipedge in TRACK. #850. (Andrea)
*    A corrector could have a length that was disregarded by TRACK. #851. (Tobias)
*    RF cavity defined with HARMON and zero LENGTH gets FREQ set to zero. 132 #858. (Tobias)
*    Prevent infinite loop when aperture is exactly on closed orbit. #874. (Tobias)
*    Fix of the thin combined function magnets in TWISS. #898. (Malte, Riccardo, Tobias)
*    Fixes in the aperture algorithm as well as the addition of dispersion. #890. (Tobias, Thys)
*    Corrected the scaling of the TWISS functions with PT in PTC. #903. (Tobias, Piotr)
*    BV flag fixed for the element CRABCAVITY to have correct behavior. #919. (Tobias, Riccardo) 

MAD-X release 5.05.02 (2019.07.26)

*    Bug fix in MAKETHIN. In version 5.05.00 and 5.05.01 correctors were sliced but the strength was not changed. #801. (Helmut, Tobias)
*    Reverted equations for PT in module EMIT -> Same partition number as before 5.05.00. #784. (Andrea)
*    x,y - rotations are now exported to SixTrack . #780. (Tobias)
*    Added several mask files to the test suite. #789. (Tobias)
*    Fixes the problem that several survey changes the comments. #788. (Tobias) 

MAD-X release 5.05.01 (2019.06.07)

*    Fixed bug that converted K3 values to K3S in MAKETHIN for octupoles. This error was introduced in the 5.05.00 release. #781. (Laurent, Tobias, Ewen) 

MAD-X release 5.05.00 (2019.05.10)
Backward incompatibilties with previous releases

*    Changes the sign of the y-rotation in the survey. #729. (Tobias, Laurent)
*    Prevents energy spread always beeing recalculated from the longitudinal emittance. #755. (Tobias) 

Bug Fixes

*    Bug fix in RF-cavity lag adjustment for no_cavity_totalpath=true. #759. (Laurent)
*    Touchek lifetime bug when negative charge was used. #754. (Tobias)
*    Fixed bug in PTC_TWISS and PTC_NORMAL that did not add deltap to the user orbit. #751. (Piotr, Laurent)
*    Prevents conversion of an aperture in case no size has been given. #736. (Tobias)
*    Bugfix for coasting beam in the IBS module. #735. (Angela)
*    Frequency was not output in SXFwrite. #734. (Tobias)
*    Fixed a problem with readtable/write tables. #720. (Laurent)
*    MKTHIN was previously disregarding disregards E1 and E2 of RBEND. #710. (Helmut, Tobias)
*    The thick quadrupole was tilted twice in TRACK when k1s was not equal to 0. #706. (Laurent)
*    Fix selecting by class=variable and class=sequence. #696. (Thomas)
*    A bug when cycling the lattice that created a negative drift with exactly the length of the sequence. #683. (Tobias) 

Added features/checks/functionalities

*    New PTC version (2019.04.18). #760. (Piotr)
*    Mapdump added to PTC. #759. (Laurent)
*    New regression tests for all elements in PTC for many configurations. #759. (Laurent)
*    Comments implemented in survey. #755. (Tobias)
*    New flag in SixTrack converter for more info in file fc.2. #749. (Veronika)
*    Fixed synchrotron radiation equations to be fully relativistic, in TWISS, TRACK, and EMIT. #762. (Andrea)
*    Fixed momentum damping equations in TWISS, TRACK, and EMIT. #762. (Andrea)
*    Adds element solenoid to module EMIT and added synrad to solenoid in TWISS. #748. (Andrea)
*    Prevents conversion of an aperture in case no size has been given. #736. (Tobias)
*    SixTrack converter: Long names for aper file. #746. (Veronika)
*    New version of MAKETHIN. #753. (Helmut, Tobias)
*    Add `chdir` command to change directory. #741. (Thomas)
*    Tilted solenoid in TWISS. #737. (Tobias, Helmut)
*    Option to update element based on changes from a parent element. #733. (Tobias)
*    Translate command does now work as a virtual drift. #727. (Tobias)
*    Added the support to retrive the errors as an internal table. #725. (Tobias)
*    Added conversion of the generalized RF-multipoles. #714. (Tobias)
*    Added support for thick RBEND in TRACK. #697. (Andrea)
*    Move command now keeps variables when the sequence is saved. #692. (Tobias)
*    Added beam sanity checks such as checking that the energy is > 0. #689. (Kyrre) 

MAD-X release 5.04.02 (2018.10.03)

*    Added a X11=yes/no switch to the build systems do disable usage of X11 on linux/mac #652. (Thomas)
*    Fix a bug in the export of the trombone matrix to Sixtrack. #656. (Joel, Tobias)
*    Increased numerical stability in PTC when TOTALTIME=true and many turns. #657. (Piotr)
*    Added full functionality to find wavelength to subtract in TOTALPATH=true mode. #663. (Piotr)
*    Fix make thin for k0=0 when angle is also defined. #666. (Tobias, Helmut)
*    Prevents the aperture model running with the wrong TWISS settings. #667. (Tobias)
*    Fix bug in racetrack aperture check. #673. (Laurent)
*    Fixed memory corruption after several TWISS commands when comments were selected. #677. (Tobias)
*    Fixed plotting when there are only solenoids in the lattice. #679. (Tobias) 

MAD-X release 5.04.01 (2018.07.10)

*    Increased the possible length of a name in the conversion to SixTrack. #620. (Tobias, James)
*    Added the possibility comments to elements. #619. (Laurent, Tobias)
*    Twiss output contains information about the beam-beam element. #622. (Tobias)
*    The BV flag is included in the twiss output. #367. (Tobias)
*    Fixed the documentation about E1 and E2 in Bending Magnet section. #618. (Laurent)
*    Corrected factorials for normal form results. #617. (Piotr)
*    RDT oscillations in ptc_twiss. #611. (Piotr)
*    Disabling stochastic for closed orbit search also in ptc_twiss and ptc_normal. #611. (Piotr)
*    Aper_offset and mech_sep now in slices after makethin. #615. (Riccardo)
*    Fix incorrect field errors in tmbend with INTERPOLATE. The full length was not used. #610. (Thomas)
*    Fix parsing negated logical after string parameter. #601. (Thomas)
*    The output to Sixtrack of the Trombone element is now in sigma and p_sigma. #605. (Tobias)
*    Synchrotron Radiation generator added. #608,#609. (Helmut, Laurent)
*    Possible to set `sectorfile=""` to disable file output. #606. (Thomas)
*    SROTATION in SURVEY is changed back to the same defintion as in 5.03.07. #642. (Thomas) 

MAD-X release 5.04.00 (2018.03.02)
Backward incompatibilties with previous releases

*    RBARC attribute of OPTION is now considered by SURVEY.
*    SROTATION angle sign changed in SURVEY for consistency with MAD-X and MAD-X PTC.
*    YROTATION (and XROTATION) properly implemented and made consistent with MAD-X PTC.
*    ELSEPARATOR kicks were swapped in H-V planes by TWISS. 

Bug Fixes

*    Fix lost support for MVAR in match #543. (Tobias)
*    Fix problem with ESAVE #501. (Thomas)
*    Fix consistency in RF-multipole kick phase when exported to SixTrack #505. (Tobias)
*    Fix ELSEPARATOR. Before a horizontal field gave a vertical kick in twiss #509. (Tobias)
*    Sign changed in SURVEY for SROTATION. #569. (Tobias)
*    Fix the replace command that changed type of other elements. #599. (Tobias)
*    Fix problem to read turn number and particle ID from TRACK table #533. (Tobias)
*    Fix type and length of aperture elements name in sixtrack export #552. (Veronika, James)
*    Fix the difference of output for VALUE command depending on compiler. #531. (Piotr)
*    Ptc_track with element_by_element=false, before it was crashing in some cases (Piotr)
*    Sign of rfcavity LAG was swapped in PTC #539. (Piotr)
*    Sign of T of closed orbit in summary table was swapped #539. (Piotr)
*    Removed absolute limit for T coordinate in ptc_tracking #548. (Piotr)
*    Fixed ptc_align #545. (Piotr)
*    6D closed orbit search with TOTALPATH=true #551. (Piotr)
*    Disabling STOCHASTIC switch during closed orbit search (was always failing) #565. (Piotr)
*    Fixed the calculation of energy loss due to quantum synchrotron radiation emission. (Andrea)
*    Fix aperture plot #490. (Thomas)
*    Fix gub in SAVE selection #492. (Thomas)
*    Fixing the unified atomic mass #518. (Laurent, Thomas)
*    Fix plotting of two points with the same x-value. #576. (Tobias)
*    Fix SURVEY to be aware of the RBARC option #564. (Laurent)
*    Fixing the last digits of the physical constants from PDG in the code. (Laurent) 

Added features/checks/functionalities

*    Thick 6d SOLENOID tracking implemented in TRACK #495 #564. (Andrea)
*    The TRANSLATION element implemented #515. (Tobias)
*    SECTORPURE option implemented for sectormaps #518. (Tobias)
*    Thick attribute for SOLENOID implemented in MAKETHIN #547. (Helmut)
*    XROTATION and YROTATION implemented. #569. (Laurent, Tobias)
*    Field errors added to thick SBENDs and QUADRUPOLEs in TRACK #583. (Andrea)
*    MULTIPOLEs and SBENDs support k0l != angle in SURVEY, TRACK, TWISS and EMIT #564. (Laurent)
*    Implemented new thick map for combined-function SBEND (and QUADRUPOLE) with h, k0, and k1 handled independently #471. (Andrea)
*    Name mangling in SixTrack export to uniquely compress names longer than 16 characters. (Andrea)
*    PTC_CREATE_LAYOUT now stops if a not supported aperture is found in any element. (Piotr)
*    Added AC dipole in PTC_TRACK and PTC_TRACKLINE (Piotr) #531. (Piotr)
*    Add an initial guess attribute for closed orbit search in PTC_TRACK #567. (Piotr)
*    Implemented roll errors of thin dipoles (EALIGN in thin lattice are very approximative) #549. (Piotr)
*    Orbit_-cT renamed to orbit_t to allow access and use within expressions #548. (Piotr)
*    PTC_SETSWITCH modifies only the parameters explicitly defined by the user (no more default update) #548. (Piotr)
*    Added seed option in PTC_SETSWITCH to set random number generator seed of PTC #548. (Piotr)
*    Added option TRACKRDTS in PTC_TWISS to produce table TWISSRDT with Resonance Driving Terms (order >=3) #548. (Piotr)
*    Automatic determination of SECTOR_NMUL parameter, speeds up calculations with exact Hamiltonian (EXACT=true) by factor 3 #559. (Piotr)
*    Updated PTC version (fix problems with chromaticities of 4th order and above) #537. (Piotr)
*    Review of coupling calculation in MAD-X, slides part 1 and part 2. (Laurent)
*    Review of SOLENOID maps and YROTATION maps for tilted solenoid-antisolenoid of FCC-ee, slides. (Laurent, Tobias)
*    Switch to C++11 (Laurent, Piotr, Thomas) 

MAD-X release 5.03.07 (2017.10.20)

*    Fix signs for Sixtrack export of Beam-Beam #467. (Riccardo)
*    Fix curvature in thick sbend #483. (Andrea)
*    Fix segfault in Best random number generators #476. (Laurent)
*    Fix in limitation of eigenvalues calculation #455. (Tobias)
*    Fix PTC_twiss quitting without tiding up memory #462. (Piotr)
*    Fix in sigma matrix initialisation in presence of coupling #425. (Tobias)
*    Fix in table read/write with duplicated lines occurring #489. (Tobias)
*    Fix problem with plotting columns from a read Twiss table #488. (Tobias)
*    Fix coupling checks #484. (Tobias, Laurent)
*    Fix documentation for time definition #470. (Tobias)
*    Fix BBOrbit flag #466. (Laurent)
*    Tests added for checking coupling in challenging setup. (Tobias)
*    Abort PTC_create_layout when matrix element encountered (was ignored). (Piotr)
*    Fix radiation in PTC_twiss and implement adequate tests on request of FCC-ee users. (Piotr)
*    Fix energy gain calculation for linac acceleration for time=true. (Piotr)
*    Add PTC_track radiation argument with stronger precedence than PTC_setswitch. (Piotr)
*    Add PSI angle (around S axis) rotation in ptc_eplacement. (Piotr)
*    Replace PTC_track normal form calculation by complex version for closed solutions. (Piotr)
*    Document ONETABLE option and format of TRACKONE table. (Piotr)
*    Add and document RECLOSS option in PTC_track. (Piotr)
*    Add ROOTNTUPLE option in PTC_track. (Piotr) 

MAD-X migrated on GitHub/Git  (2017.06.23)

*    MAD-X new repository on GitHub. 

MAD-X release 5.03.06 (2017.05.29)

*    Bug fix in dipedge element map [Laurent]. 

MAD-X release 5.03.05 (2017.05.17)

*    Bug fix in sigma matrix initialization [Irina]. 

MAD-X release 5.03.04 (2017.05.15)

*    Radiation bux fixes (including quantum), now working properly [Andrea].
*    Matching bux fixes, now working with the chrom option (ticket #414) [Laurent].
*    Beam attribute access bufg fix (nasty side effect, ticket #422) [Laurent].
*    New random number generator with support for 10 streams [Laurent].
*    New track quantum seed option added [Laurent]. 

MAD-X release 5.03.00 (2017.04.14)

*    Adding setvars_const, setvars_knob, fill_knob commands and fill->scale parameter [Riccardo].
*    Computation of higher order alpha_c with ptc_twiss, time=true, icase=56 [Piotr].
*    Fix of issue in ptc_twiss big slow down (ticket #417) [Piotr].
*    Scaling of dispersions and chromaticities with beta for ptc_twiss with time=true [Piotr].
*    Cleaner higher order compaction factor calculation [Piotr].
*    Bug fix bug (ticket #402) where T and PT were not checked if within user specified box [Piotr].
*    Bug fix to get Hamiltonian terms from ptc_normal with the new PTC [Piotr]. 

MAD-X new release numbering  (2017.04.14)

*    MAD-X major number represents (pro) release number, minor number represents #bug fixes. 

MAD-X release 5.02.13 (2016.12.20)

*    Setvars_knob and fill_knob new commands plus testsuite [Riccardo]
*    Documentation updates for tickets overcome and answers [Irina, Andrea, Laurent]
*    KSL0 error added to kickers implementation (ticket #339) [Irina, Andrea]
*    Track maxaper limits checks enabled and extended to longitudinal [Irina]
*    Update of cosmux and cosmuy checks (now 3 methods are implemented for xchecks) [Irina]
*    Twiss chrom bug fix in case of failure with first deltap trial (wrong chromaticity) [Irina]
*    Inaccurate chromatic functions calculation for coupled lattice, warning added in Twiss [Irina]
*    Update of tests to follow bugs fixes and cleanup [Irina, Laurent]
*    Bug fix in factorial for anharmonicities in ptc_normal [Piotr]
*    Bug fix (ticket #394) in TFS table for error (unchecked column access) [Piotr]
*    K1 implemented in thick sbend in track [Andrea]
*    Makethin does not share anymore last slice of thick quad to allow proper name (ticket #329) [Andrea]
*    Bug fix and gards added in setvars/fill commands familly [Laurent]
*    Bug fix for angles and bv flag in Survey command, update of TI8 sequence [Laurent]
*    Update formula for multipole angle and tilt calculation, update documentation [Laurent]
*    Survey command updated for pitfalls in documentation [Laurent]
*    Lrad and tilt added to monitors (not used, ticket #361) [Laurent]
*    Bug fix about invalid use of TKicker by Track and Twiss (ticket #398) [Laurent]
*    Space charge final update in Track and Twiss [Frank, Laurent]
*    Bug fix in dipedge tilt in Track and translation to PTC [Frank, Laurent]
*    Bug fix in beam parameter, relativistic beta missing in one update case [Frank, Laurent]
*    Ongoing integration of general multipole (was Combined Function Magnet) TBC [Malte, Laurent]
*    libmadx-linux32 not built anymore due dependencies to missing libX11.so 32 bit [Laurent] 

MAD-X release 5.02.12 (2016.10.15)

*    Bug fix for wrong logic in coupling [Irina]
*    Bug fix in Track command [Laurent]
*    Bug fix in parser [Laurent]
*    Keeptrack option for Run command [Laurent] 

MAD-X release 5.02.11 (2016.10.11)

*    Coupling fixes (see slides) [Irina]
*    Secondary method to cross-check fractional part of the tunes (see slides) [Irina]
*    Sigma matrix in twiss (i.e. sig11..sig66) [Irina]
*    Correction of a logic error in max_mult_ord settings in c6t [Laurent]
*    New PTC and PTC tests update [Piotr]
*    Bug fixes caused by all variables initialized within declaration in Fortran 95 (#389) [Laurent,Irina,Andrea]
*    Committing F. Schmidt's Space-Charge update [Andrea]
*    Minor improvement of the build system [Laurent]
*    Improvement of the CMake build system [Yngve]
*    Add flag to export all markers to Sixtrack (#391) [Laurent]
*    Some bug fixes (e.g. #388, #394) [Laurent, Piotr] 

MAD-X release 5.02.10 (2016.04.23)

*    New Twiss command attribute 'sectoracc' to twiss command for saving composed maps in sectormap [Laurent]
*    Bug fix in coupling for phase advance computation [Laurent, Irina]
*    Default beam 'et', 'sigt', 'sige' set in dictionnary instead of beam update [Laurent]
*    New fix point for probe initialization (affected if RF and/or dp) (#385) [Laurent]
*    Cleaning of modules initialization: twiss, track, emit, touschek, dynap, plot, ibs [Laurent]
*    New 'tabindex' feature to find the index of a string in a table column [Laurent]
*    More attributes copied by makethin to generated thin sequence [Helmut] 

MAD-X release 5.02.09 (2016.04.06)

*    Many bug fixes.
*    Bug fix for one turn map initialization in many modules (see slides) [Laurent]
*    Bug fix for RF initialization [Laurent]
*    Bug fix in coupling [Irina]
*    Bug fix in micado [Ghislain]
*    Flip mode implementation for strong coupling [Irina]
*    New scheme to synchronize C vs Fortran I/O [Laurent]
*    New Windows (10) build using msys2+mingw32+mingw64 platforms [Laurent]
*    New Linux (Ubuntu 14.04) build for standalone linux [Laurent]
*    Garbage collector natively build and tested on Windows [Laurent]
*    Default release set to GNU 64 bit [Laurent]
*    ... 

MAD-X new Linux and Windows platforms  (2016.04.06)

*    MAD-X is now built and tested on Linux Ubuntu and Windows 10 (+msys2) VMs. 

MAD-X release 5.02.08 (2015.11.23)

*    Release for stabilization before pro-release.
*    Policy changes for release delivery.
*    LHC and HLHC user-cases added to the test system.
*    Backward compatibility checks with LHC/HLLHC studies (see slides).
*    Some bug fixes and new feature (see tickets).
*    Makethin thick slices improved, default teapot and makedipedge [Helmut]
*    Documentation correction and improvement [Ghislain].
*    ... 

MAD-X release 5.02.07 (2015.09.25)

*    Many bug fixes.
*    Makethin thick slices improved, default teapot [Helmut, Laurent]
*    Documentation correction and improvement [Ghislain].
*    ... 

MAD-X release 5.02.06 (2015.06.26)

*    Many bug fixes.
*    New makethin [Helmut, Laurent]
*    Documentation correction and improvement [Ghislain].
*    ... 

MAD-X release 5.02.05 (2015.04.10)

*    Minor release with some bug fix.
*    Documentation correction and improvement [Ghislain].
*    ... 

MAD-X release 5.02.04 (2014.11.14)

*    Minor release with some bug fix.
*    Documentation with some correction, including figure 24.3 [Ghislain].
*    Portable commanded copyfile, renamefile, removefile added and unified [Laurent].
*    Column "Keyword" added to survey table for consistency with twiss table (#308) [Laurent]. 

MAD-X new users guide  (2014.10.16)

*    MAD-X new user's guide written in LaTeX and distributed in PDF. 

MAD-X release 5.02.03 (2014.10.15)

*    New Latex/PDF MAD-X userguide (draft) [Ghislain].
*    Add support for sanitizer (#305) [Laurent].
*    Add support for openmp (#299) [Laurent].
*    Bug fix in element management (#296) [Laurent, Piotr].
*    Bug fix in garbage collector with Intel (#307) [Laurent].
*    Bug fix in survey to process skew component in multipole (#306) [Laurent].
*    Bug fix in the initial orbit input for twiss (#271) [Ghislain].
*    Bug fix in use_macro on MacOSX 10.9 (#298) [Piotr].
*    Bug fix in ptc_twiss to handle beam-beam and space charge correctly (#295) [Piotr].
*    Bug fix to remove dummy CELL command (#303) [Laurent].
*    Bug fix in seqedit cycle after makethin (#286) [Ghislain].
*    Bug fix in LINE with multiplied elements (#290) [Ghislain, Piotr].
*    Bug fix in ptc_create_universe crashing ptc_create_layout (#287) [Piotr].
*    Bug fix and improvement in synchrotron radiation calculation (#266) [Ghislain].
*    Bug fix in ptc_normal crashing MAD-X (#282) [Laurent].
*    Error reporting improvement when PTC does not find the closed orbit (#152) [Piotr].
*    Improvement of makethin module (#283) [Helmut].
*    Many minor bugs fix. 

MAD-X release 5.02.02 (2014.07.04)

*    Many bugs fix.
*    ...
*    See also the slides of MAD-X meeting on indico. 

MAD-X memory management  (2014.07.04)

*    MAD-X memory footprint is halved for standard (HL-LHC) studies. 

MAD-X release 5.02.01 (2014.04.25)

*    Minor bugs fix.
*    Bug fix in plot (#261) [Laurent, Ghislain].
*    Bug fix in table ranges (#262) [Laurent].
*    Bug fix in numdiff [Laurent].
*    Bug fix in refpos (#257) [Ghislain].
*    Bug fix in circular sequence [Ghislain].
*    Bug fix in threader (#260) [Ghislain].
*    Bug fix in beam update (#259) [Ghislain]. 

MAD-X release 5.02.00 (2014.03.05)

*    2nd production release since 2011.06.
*    Minor bugs fix.
*    Major bug fix in aperture when loading offset (ticket #254) [Laurent, Ghilain]
*    New printf-like command [Laurent]. 

MAD-X pro-release  (2014.03.05)

*    MAD-X pro-release 5.02.00 becomes the default on lxplus. 

MAD-X 64 bit  (2014.02.13)

*    MAD-X 64 bit becomes the default on lxplus (see releases). 

MAD-X release 5.01.06 (2014.01.31)

*    Minor bugs fix.
*    Summary of orbit correction changed (bug fix) [Ghislain].
*    ptc_twiss: [Piotr]
*    implementation of map_table in ptc_twiss the same way as in ptc_normal.
*    support for initial condition specification (map table).
*    implementation of min and max for betas and dispersion, and min, max and rms for the orbit. 

MAD-X release 5.01.05 (2013.12.20)

*    Minor bugs fix.
*    Upgrade to last PTC release [Laurent].
*    Random crash corrected (thanks to new win-gnu arch) [Laurent]. 

MAD-X release 5.01.04 (2013.12.13)

*    Minor bugs fix.
*    Nightly tests and build fully operational [Laurent].
*    Cross platforms build and test reports completed [Laurent].
*    Update of the physical constants published by the PDG in 2012 [Ghislain].
*    New support for platform Gnu windows (11 platforms supported now) [Laurent].
*    Fixed numerical instability in Touschek module and improved parameter handling for integration [Ghislain]. 

MAD-X release 5.01.03 (2013.11.05)

*    Many bugs fix.
*    Nested sequence bug corrected [Ghislain].
*    New thick quadrupole element for tracking [Andrea].
*    New makethin module (in C++) [Helmut, Laurent].
*    New static space charge code in track (expert only) [Frank, Laurent].
*    MAD-X on indico [Laurent]. 

MAD-X nightly builds and tests  (2013.10.21)

*    MAD-X moved to a dedicated server for builds and tests with nightly reports (see development). 

MAD-X on Indico  (2013.10.07)

*    MAD-X events are now archived on Indico (see events). 

MAD-X numdiff manual  (2013.08.14)

*    MAD-X numdiff tool has been extended for data validation and the manual is published. 

MAD-X release 5.01.02 (2013.05.30)

*    Minor bugs fix.
*    Numdiff manual added [Laurent].
*    New multi-pass tracking feature [Frank, Laurent]. 

MAD-X release 5.01.01 (2013.05.10)

*    Minor bugs fix.
*    Garbage collector to fix memory leaks [Laurent].
*    Full 64 bit support (better than 32 bit now) [Laurent]. 

MAD-X memory management  (2013.05.10)

*    MAD-X dynamic memory allocations are now garbage collected (see roadmap). 

MAD-X release 5.01.00 (2013.03.13)

*    1st production release since 2011.06.
*    Minor bugs fix.
*    New numdiff [Laurent].
*    Input script as argument (input redirection no more needed) [Laurent]. 

MAD-X pro-release  (2013.03.13)

*    MAD-X pro-release 5.01.00 becomes the default on lxplus. 

MAD-X release 5.00.20 (2013.02.20)

*    Many bugs fix.
*    New table API (major bug fix) [Laurent].
*    Upgrade to last PTC release [Laurent].
*    New removefile, renamefile, appendfile commands for script portability [Laurent].
*    Tests system completed for Windows, Linux and MacOSX in 32 and 64 bit (234 tests) [Laurent, Andrea, Piotr]. 

MAD-X on OS X 10.5  (2013.01.08)

*    MAD-X users running OS X 10.5 (Leopard) must use madx releases 5.00.19 or later. 

MAD-X project  (2012.11.20)

*    The MAD-X status and future plans slides (restricted access). 

MAD-X module keepers  (2012.10.22)

*    MAD-X module keepers team as been updated (see contributors). 

MAD-X test suite  (2012.09.26)

*    MAD-X first test suite now contains 37 tests (see roadmap). 

MAD-X at ICAP'12  (2012.08.23)

*    The MAD-X presentation and paper presented at the ICAP'12 conference. 

MAD-X shared library  (2012.06.27)

*    MAD-X is now provided as a standalone library for scripting languages (see roadmap). 

MAD-X test system  (2012.06.08)

*    MAD-X is now equiped with a new test system (see roadmap). 

MAD-X numdiff  (2012.05.30)

*    MAD-X provides now the numdiff tool to compare numeric files (see roadmap). 

MAD-X/PTC library  (2012.03.14)

*    MAD-X provides now PTC as a standalone library for lxplus users (see roadmap). 

MAD-X workspace  (2012.03.07)

*    MAD-X has now its own workspace on AFS (see roadmap). 

MAD-X examples  (2012.03.05)

*    The examples has been moved to SVN (see roadmap). 

MAD-X trac system  (2012.02.01)

*    The trac system is now fully operational (see roadmap). 

MAD-X user's guide  (2012.01.18)

*    The user's guide has been moved to SVN (see roadmap). 

MAD-X website  (2011.12.08)

*    The new website (this website) has been published (see roadmap). 

MAD-X build system  (2011.11.28)

*    MAD-X is now equiped with a brand new build system (see roadmap). 

First 10 years of PTC  (2011.11.09)

*    The PTC 1/2 day workshop (restricted access). 

MAD-X 1st contact  (2011.05.31)

*    The MAD-X project overview presented at the BE-ABP-LCU meeting (restricted access). 




