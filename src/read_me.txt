2025-02-18
There are too many missing values in glwd.nc, especially in the study region. On a second thought, I replaced the m_fwlnd values with vagetation-specific parameters.
This can be easily changed (back) in data_buffer_mod.f90. Notice: I added a parameter fwlnd in the subroutine which can be removed if not needed.
Updated carbon_cycle_mod and bubble_mod according to ALBM 2024.Nov update on Github'
2025-02-19
According to Tan's suggestion and comments, SOCref is back to the value of 21.98 kgC m-2 while m_soc uses the vegetation-specific value instead of reading from HWSD v1.2.
The changes including bg_const_mod (adding 0218 changes into this mod instead) and data_buffer_mod.

2025-02-19
Decide to incorporate Tan et al., 2024
(1) CH4 oxidation. Change Ko2 and Kch4 to BetaCH4 and LamO2 in shr_param_mod and bg_utilities_mod,
Changes the value in opt_BLH1 for test.

2025-04-02
Find a bug in diagenesis.mod.

2025-04-14
change the calculation for CH4 spread from sediment to water to the original way.
Change the methane generating rate from 2.5*1d-9 to 1d-9 in the function(undo)
Change the function for methane oxidation

2025-04-15
change the function for methane generation and oxidation

2025-04-16
Change bubble content ratio
Change back winter bubble raio
change oq10 from 1.7 to 3.5

2025-04-20
Change AOM ratio from 1.5 to 2.7
Change back Ebch4 to a variable to depth

2025-04-22
Fix top=watertopindex instead of 1, though they are not yedoma.
Change back oq10

2025-04-23
        Change bubble flux(the last subroutine)

2025-04-24
Change back of 20250422
Change Photosynthesis: ftemp=max(e8,..) instead of max(0.0,..)
AOMQ10&PQ10 changed, all according to Xu(back)
Change anaerobic oxidation occurance condition and the ratio to production rate
Change aerobic oxidation function (1d-4 to 1d-3)

2025-09-03
(1) diagenesis_mod: Shift labile from 0.05 to 1
(2) bg_utilities_mod: Kill phBeta and all parts inc. phBeta
(3) shr_param_mod: Re-numbering; only keep those effective. Add "Kan"(0.15~0.27) as sensitivity parameter.
(4) bg_utilities_mod: Add "Kan" instead of fixed 0.15 when calculating anaerobic oxidation.
(5) m_Qgw: Use month mean of pr instead of daily. Add gw_file in namelist.
(6) data_buffer_mod: Read gw_file.
(7) shr_ctrl_mod: Add gw_file.
(8) carbon_cycle_mod: Introduce in "isDayNode" in CalcCarbonCycleRates(), thus CarbonModuleSetup(), as well as in sim_coupler_mod. So that the dark fixation is only calculated once a day.
(9) diagenesis_mod: Change back in subroutine CalcDiagenesisRates(), about those prtot and Ebch4().
(10) read_data_mod: Add gw_file into run_data group.

2025-09-04
(1) diagenesis_mod: Include capillary effects. Ebullition is triggered shen the sum of dissolved-gas partial pressures exceeds sum(atmospheric+hydrostatic+capillary).
(2) diagenesis_mod: Fix the proprtion of different gases in the total pressure and bubbles, including N2, O2, CO2 and CH4.

