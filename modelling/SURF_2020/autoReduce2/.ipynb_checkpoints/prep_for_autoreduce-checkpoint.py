# imports
from libsbml import *
import sys
import numpy as np
from auto_reduce import *
from auto_reduce.utils import get_ODE
import matplotlib.pyplot as plt
from auto_reduce.utils import reduce

# set num species and params
n = 68
num_params = 80

r0_k = 15.0
r0_1_k = 1.5
r2_k = 36000.0
r3_k = 15.0
r3_1_k = 1.5
r5_k = 15.0
r5_1_k = 1.5
r7_k = 36000.0
r8_k = 15.0
r8_1_k = 1.5
r10_k = 15.0
r10_1_k = 1.5
r12_k = 36000.0
r13_k = 15.0
r13_1_k = 1.5
r15_k = 15.0
r15_1_k = 1.5
r17_k = 36000.0
r18_k = 15.0
r18_1_k = 1.5
r20_k = 15.0
r20_1_k = 1.5
r22_k = 36000.0
r23_k = 15.0
r23_1_k = 1.5
r25_k = 15.0
r25_1_k = 1.5
r27_k = 36000.0
r28_k = 15.0
r28_1_k = 1.5
r30_k = 15.0
r30_1_k = 1.5
r32_k = 36000.0
r33_k = 15.0
r33_1_k = 1.5
r35_k = 15.0
r35_1_k = 1.5
r37_k = 36000.0
r38_k = 15.0
r38_1_k = 1.5
r40_k = 15.0
r40_1_k = 1.5
r42_k = 36000.0
r43_k = 15.0
r43_1_k = 1.5
r45_k = 15.0
r45_1_k = 1.5
r47_k = 36000.0
r48_k = 15.0
r48_1_k = 1.5
r50_k = 15.0
r50_1_k = 1.5
r52_k = 36000.0
r53_k = 15.0
r53_1_k = 1.5
r55_k = 15.0
r55_1_k = 1.5
r57_k = 36000.0
r58_k = 15.0
r58_1_k = 1.5
r60_k = 15.0
r60_1_k = 1.5
r62_k = 36000.0
r63_k = 15.0
r63_1_k = 1.5
r65_k = 15.0
r65_1_k = 1.5
r67_k = 36000.0
r68_k = 15.0
r68_1_k = 1.5
r70_k = 15.0
r70_1_k = 1.5
r72_k = 36000.0
r73_k = 15.0
r73_1_k = 1.5
r75_k = 15.0
r75_1_k = 1.5
r77_k = 2.0
r78_k = 15.0
r78_1_k = 1.5


P = np.zeros(num_params)
# Now set all parameters
i = 0
P[i] = r0_k
i += 1
P[i] = r0_1_k = 1.5
i += 1
P[i] = r2_k = 36000.0
i += 1
P[i] = r3_k = 15.0
i += 1
P[i] = r3_1_k = 1.5
i += 1
P[i] = r5_k = 15.0
i += 1
P[i] = r5_1_k = 1.5
i += 1
P[i] = r7_k = 36000.0
i += 1
P[i] = r8_k = 15.0
i += 1
P[i] = r8_1_k = 1.5
i += 1
P[i] = r10_k = 15.0
i += 1
P[i] = r10_1_k = 1.5
i += 1
P[i] = r12_k = 36000.0
i += 1
P[i] = r13_k = 15.0
i += 1
P[i] = r13_1_k = 1.5
i += 1
P[i] = r15_k = 15.0
i += 1
P[i] = r15_1_k = 1.5
i += 1
P[i] = r17_k = 36000.0
i += 1
P[i] = r18_k = 15.0
i += 1
P[i] = r18_1_k = 1.5
i += 1
P[i] = r20_k = 15.0
i += 1
P[i] = r20_1_k = 1.5
i += 1
P[i] = r22_k = 36000.0
i += 1
P[i] = r23_k = 15.0
i += 1
P[i] = r23_1_k = 1.5
i += 1
P[i] = r25_k = 15.0
i += 1
P[i] = r25_1_k = 1.5
i += 1
P[i] = r27_k = 36000.0
i += 1
P[i] = r28_k = 15.0
i += 1
P[i] = r28_1_k = 1.5
i += 1
P[i] = r30_k = 15.0
i += 1
P[i] = r30_1_k = 1.5
i += 1
P[i] = r32_k = 36000.0
i += 1
P[i] = r33_k = 15.0
i += 1
P[i] = r33_1_k = 1.5
i += 1
P[i] = r35_k = 15.0
i += 1
P[i] = r35_1_k = 1.5
i += 1
P[i] = r37_k = 36000.0
i += 1
P[i] = r38_k = 15.0
i += 1
P[i] = r38_1_k = 1.5
i += 1
P[i] = r40_k = 15.0
i += 1
P[i] = r40_1_k = 1.5
i += 1
P[i] = r42_k = 36000.0
i += 1
P[i] = r43_k = 15.0
i += 1
P[i] = r43_1_k = 1.5
i += 1
P[i] = r45_k = 15.0
i += 1
P[i] = r45_1_k = 1.5
i += 1
P[i] = r47_k = 36000.0
i += 1
P[i] = r48_k = 15.0
i += 1
P[i] = r48_1_k = 1.5
i += 1
P[i] = r50_k = 15.0
i += 1
P[i] = r50_1_k = 1.5
i += 1
P[i] = r52_k = 36000.0
i += 1
P[i] = r53_k = 15.0
i += 1
P[i] = r53_1_k = 1.5
i += 1
P[i] = r55_k = 15.0
i += 1
P[i] = r55_1_k = 1.5
i += 1
P[i] = r57_k = 36000.0
i += 1
P[i] = r58_k = 15.0
i += 1
P[i] = r58_1_k = 1.5
i += 1
P[i] = r60_k = 15.0
i += 1
P[i] = r60_1_k = 1.5
i += 1
P[i] = r62_k = 36000.0
i += 1
P[i] = r63_k = 15.0
i += 1
P[i] = r63_1_k = 1.5
i += 1
P[i] = r65_k = 15.0
i += 1
P[i] = r65_1_k = 1.5
i += 1
P[i] = r67_k = 36000.0
i += 1
P[i] = r68_k = 15.0
i += 1
P[i] = r68_1_k = 1.5
i += 1
P[i] = r70_k = 15.0
i += 1
P[i] = r70_1_k = 1.5
i += 1
P[i] = r72_k = 36000.0
i += 1
P[i] = r73_k = 15.0
i += 1
P[i] = r73_1_k = 1.5
i += 1
P[i] = r75_k = 15.0
i += 1
P[i] = r75_1_k = 1.5
i += 1
P[i] = r77_k = 2.0
i += 1
P[i] = r78_k = 15.0
i += 1
P[i] = r78_1_k = 1.5
i += 1

params_values = P.copy()
timepoints_ode = np.linspace(0,75,100)

error_tol = 100
nstates_tol = 1


# Set initial conditions
# want all enzymes = 0.15, glucose, atp, nadp, maybe pi = 30, atpase_init = 2.0
x_init = np.zeros(n)

# set up init for all enzymes
enz_list = [0, 7, 11, 15, 19, 25, 30, 33, 37, 41, 45, 49, 53, 57, 61]
# indices for gluc, atp, nadp, pi
mol_list = [2, 1, 20, 26]

for i,item in enumerate(x_init):
    if i in enz_list:
        x_init[i] = 0.15
        
    # atpase init cond
    if i == 65:
        x_init[i] = 2.
        
    # set up other molecules
    if i in mol_list:
        x_init[i] = 30


x,f,P = system.load_ODE_model(n, num_params)

params = P

enzyme_hex = x[0]
metabolite_atp = x[1]
molecule_glucose = x[2]
molecule_g6p = x[3]
metabolite_adp = x[4]
complex_enzyme_hex_metabolite_atp_molecule_glucose = x[5]
complex_enzyme_hex_metabolite_adp_molecule_g6p = x[6]
enzyme_pgi = x[7]
molecule_f6p = x[8]
complex_enzyme_pgi_molecule_g6p = x[9]
complex_enzyme_pgi_molecule_f6p = x[10]
enzyme_pfk = x[11]
molecule_f16p = x[12]
complex_enzyme_pfk_metabolite_atp_molecule_f6p = x[13]
complex_enzyme_pfk_metabolite_adp_molecule_f16p = x[14]
enzyme_ald_tpi = x[15]
molecule_g3p = x[16]
complex_enzyme_ald_tpi_molecule_f16p = x[17]
complex_enzyme_ald_tpi_2x_molecule_g3p = x[18]
enzyme_gapN = x[19]
metabolite_nadp = x[20]
molecule_3pg = x[21]
metabolite_nadph = x[22]
complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p = x[23]
complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg = x[24]
enzyme_mGapDH = x[25]
metabolite_pi = x[26]
molecule_13bpg = x[27]
complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p = x[28]
complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg = x[29]
enzyme_pgk = x[30]
complex_enzyme_pgk_metabolite_adp_molecule_13bpg = x[31]
complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg = x[32]
enzyme_pgm = x[33]
molecule_2pg = x[34]
complex_enzyme_pgm_2x_molecule_3pg = x[35]
complex_enzyme_pgm_2x_molecule_2pg = x[36]
enzyme_eno = x[37]
molecule_pep = x[38]
complex_enzyme_eno_2x_molecule_2pg = x[39]
complex_enzyme_eno_2x_molecule_pep = x[40]
enzyme_pyk = x[41]
molecule_pyruvate = x[42]
complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep = x[43]
complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate = x[44]
enzyme_alsS = x[45]
molecule_acetolac = x[46]
complex_enzyme_alsS_2x_molecule_pyruvate = x[47]
complex_enzyme_alsS_molecule_acetolac = x[48]
enzyme_IlvC = x[49]
molecule_23dih3mebut = x[50]
complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac = x[51]
complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut = x[52]
enzyme_IlvD = x[53]
molecule_3me2oxo = x[54]
complex_enzyme_IlvD_molecule_23dih3mebut = x[55]
complex_enzyme_IlvD_molecule_3me2oxo = x[56]
enzyme_kivD = x[57]
molecule_isobutanal = x[58]
complex_enzyme_kivD_molecule_3me2oxo = x[59]
complex_enzyme_kivD_molecule_isobutanal = x[60]
enzyme_yahk = x[61]
molecule_isobutanol = x[62]
complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal = x[63]
complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol = x[64]
enzyme_atpase = x[65]
complex_enzyme_atpase_metabolite_atp = x[66]
complex_enzyme_atpase_metabolite_adp_metabolite_pi = x[67]


r0 = r0_k * enzyme_hex * metabolite_atp * molecule_glucose
r0_1 = r0_1_k * complex_enzyme_hex_metabolite_atp_molecule_glucose
r2 = r2_k * complex_enzyme_hex_metabolite_atp_molecule_glucose
r3 = r3_k * complex_enzyme_hex_metabolite_adp_molecule_g6p
r3_1 = r3_1_k * enzyme_hex * molecule_g6p * metabolite_adp
r5 = r5_k * enzyme_pgi * molecule_g6p
r5_1 = r5_1_k * complex_enzyme_pgi_molecule_g6p
r7 = r7_k * complex_enzyme_pgi_molecule_g6p
r8 = r8_k * complex_enzyme_pgi_molecule_f6p
r8_1 = r8_1_k * enzyme_pgi * molecule_f6p
r10 = r10_k * enzyme_pfk * metabolite_atp * molecule_f6p
r10_1 = r10_1_k * complex_enzyme_pfk_metabolite_atp_molecule_f6p
r12 = r12_k * complex_enzyme_pfk_metabolite_atp_molecule_f6p
r13 = r13_k * complex_enzyme_pfk_metabolite_adp_molecule_f16p
r13_1 = r13_1_k * enzyme_pfk * molecule_f16p * metabolite_adp
r15 = r15_k * enzyme_ald_tpi * molecule_f16p
r15_1 = r15_1_k * complex_enzyme_ald_tpi_molecule_f16p
r17 = r17_k * complex_enzyme_ald_tpi_molecule_f16p
r18 = r18_k * complex_enzyme_ald_tpi_2x_molecule_g3p
r18_1 = r18_1_k * enzyme_ald_tpi * pow(molecule_g3p, 2)
r20 = r20_k * enzyme_gapN * pow(metabolite_nadp, 2) * pow(molecule_g3p, 2)
r20_1 = r20_1_k * complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p
r22 = r22_k * complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p
r23 = r23_k * complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg
r23_1 = r23_1_k * enzyme_gapN * pow(molecule_3pg, 2) * pow(metabolite_nadph, 2)
r25 = r25_k * enzyme_mGapDH * metabolite_pi * pow(metabolite_nadp, 2) * pow(molecule_g3p, 2)
r25_1 = r25_1_k * complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p
r27 = r27_k * complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p
r28 = r28_k * complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg
r28_1 = r28_1_k * enzyme_mGapDH * molecule_13bpg * pow(metabolite_nadph, 2)
r30 = r30_k * enzyme_pgk * metabolite_adp * molecule_13bpg
r30_1 = r30_1_k * complex_enzyme_pgk_metabolite_adp_molecule_13bpg
r32 = r32_k * complex_enzyme_pgk_metabolite_adp_molecule_13bpg
r33 = r33_k * complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg
r33_1 = r33_1_k * enzyme_pgk * pow(molecule_3pg, 2) * metabolite_atp
r35 = r35_k * enzyme_pgm * pow(molecule_3pg, 2)
r35_1 = r35_1_k * complex_enzyme_pgm_2x_molecule_3pg
r37 = r37_k * complex_enzyme_pgm_2x_molecule_3pg
r38 = r38_k * complex_enzyme_pgm_2x_molecule_2pg
r38_1 = r38_1_k * enzyme_pgm * pow(molecule_2pg, 2)
r40 = r40_k * enzyme_eno * pow(molecule_2pg, 2)
r40_1 = r40_1_k * complex_enzyme_eno_2x_molecule_2pg
r42 = r42_k * complex_enzyme_eno_2x_molecule_2pg
r43 = r43_k * complex_enzyme_eno_2x_molecule_pep
r43_1 = r43_1_k * enzyme_eno * pow(molecule_pep, 2)
r45 = r45_k * enzyme_pyk * pow(metabolite_adp, 2) * pow(molecule_pep, 2)
r45_1 = r45_1_k * complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep
r47 = r47_k * complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep
r48 = r48_k * complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate
r48_1 = r48_1_k * enzyme_pyk * pow(molecule_pyruvate, 2) * pow(metabolite_atp, 2)
r50 = r50_k * enzyme_alsS * pow(molecule_pyruvate, 2)
r50_1 = r50_1_k * complex_enzyme_alsS_2x_molecule_pyruvate
r52 = r52_k * complex_enzyme_alsS_2x_molecule_pyruvate
r53 = r53_k * complex_enzyme_alsS_molecule_acetolac
r53_1 = r53_1_k * enzyme_alsS * molecule_acetolac
r55 = r55_k * enzyme_IlvC * metabolite_nadph * molecule_acetolac
r55_1 = r55_1_k * complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac
r57 = r57_k * complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac
r58 = r58_k * complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut
r58_1 = r58_1_k * enzyme_IlvC * molecule_23dih3mebut * metabolite_nadp
r60 = r60_k * enzyme_IlvD * molecule_23dih3mebut
r60_1 = r60_1_k * complex_enzyme_IlvD_molecule_23dih3mebut
r62 = r62_k * complex_enzyme_IlvD_molecule_23dih3mebut
r63 = r63_k * complex_enzyme_IlvD_molecule_3me2oxo
r63_1 = r63_1_k * enzyme_IlvD * molecule_3me2oxo
r65 = r65_k * enzyme_kivD * molecule_3me2oxo
r65_1 = r65_1_k * complex_enzyme_kivD_molecule_3me2oxo
r67 = r67_k * complex_enzyme_kivD_molecule_3me2oxo
r68 = r68_k * complex_enzyme_kivD_molecule_isobutanal
r68_1 = r68_1_k * enzyme_kivD * molecule_isobutanal
r70 = r70_k * enzyme_yahk * metabolite_nadph * molecule_isobutanal
r70_1 = r70_1_k * complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal
r72 = r72_k * complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal
r73 = r73_k * complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol
r73_1 = r73_1_k * enzyme_yahk * molecule_isobutanol * metabolite_nadp
r75 = r75_k * enzyme_atpase * metabolite_atp
r75_1 = r75_1_k * complex_enzyme_atpase_metabolite_atp
r77 = r77_k * complex_enzyme_atpase_metabolite_atp
r78 = r78_k * complex_enzyme_atpase_metabolite_adp_metabolite_pi
r78_1 = r78_1_k * enzyme_atpase * metabolite_adp * metabolite_pi


func_array = np.array([ + (-r0) + (r0_1) + (r3) + (-r3_1),
       + (-r0) + (r0_1) + (-r10) + (r10_1) + (r33) + (-r33_1) + ((2.0)*r48) + (-(2.0)*r48_1) + (-r75) + (r75_1),
       + (-r0) + (r0_1),
       + (r3) + (-r3_1) + (-r5) + (r5_1),
       + (r3) + (-r3_1) + (r13) + (-r13_1) + (-r30) + (r30_1) + (-(2.0)*r45) + ((2.0)*r45_1) + (r78) + (-r78_1),
       + (r0) + (-r0_1) + (-r2),
       + (r2) + (-r3) + (r3_1),
       + (-r5) + (r5_1) + (r8) + (-r8_1),
       + (r8) + (-r8_1) + (-r10) + (r10_1),
       + (r5) + (-r5_1) + (-r7),
       + (r7) + (-r8) + (r8_1),
       + (-r10) + (r10_1) + (r13) + (-r13_1),
       + (r13) + (-r13_1) + (-r15) + (r15_1),
       + (r10) + (-r10_1) + (-r12),
       + (r12) + (-r13) + (r13_1),
       + (-r15) + (r15_1) + (r18) + (-r18_1),
       + ((2.0)*r18) + (-(2.0)*r18_1) + (-(2.0)*r20) + ((2.0)*r20_1) + (-(2.0)*r25) + ((2.0)*r25_1),
       + (r15) + (-r15_1) + (-r17),
       + (r17) + (-r18) + (r18_1),
       + (-r20) + (r20_1) + (r23) + (-r23_1),
       + (-(2.0)*r20) + ((2.0)*r20_1) + (-(2.0)*r25) + ((2.0)*r25_1) + (r58) + (-r58_1) + (r73) + (-r73_1),
       + ((2.0)*r23) + (-(2.0)*r23_1) + ((2.0)*r33) + (-(2.0)*r33_1) + (-(2.0)*r35) + ((2.0)*r35_1),
       + ((2.0)*r23) + (-(2.0)*r23_1) + ((2.0)*r28) + (-(2.0)*r28_1) + (-r55) + (r55_1) + (-r70) + (r70_1),
       + (r20) + (-r20_1) + (-r22),
       + (r22) + (-r23) + (r23_1),
       + (-r25) + (r25_1) + (r28) + (-r28_1),
       + (-r25) + (r25_1) + (r78) + (-r78_1),
       + (r28) + (-r28_1) + (-r30) + (r30_1),
       + (r25) + (-r25_1) + (-r27),
       + (r27) + (-r28) + (r28_1),
       + (-r30) + (r30_1) + (r33) + (-r33_1),
       + (r30) + (-r30_1) + (-r32),
       + (r32) + (-r33) + (r33_1),
       + (-r35) + (r35_1) + (r38) + (-r38_1),
       + ((2.0)*r38) + (-(2.0)*r38_1) + (-(2.0)*r40) + ((2.0)*r40_1),
       + (r35) + (-r35_1) + (-r37),
       + (r37) + (-r38) + (r38_1),
       + (-r40) + (r40_1) + (r43) + (-r43_1),
       + ((2.0)*r43) + (-(2.0)*r43_1) + (-(2.0)*r45) + ((2.0)*r45_1),
       + (r40) + (-r40_1) + (-r42),
       + (r42) + (-r43) + (r43_1),
       + (-r45) + (r45_1) + (r48) + (-r48_1),
       + ((2.0)*r48) + (-(2.0)*r48_1) + (-(2.0)*r50) + ((2.0)*r50_1),
       + (r45) + (-r45_1) + (-r47),
       + (r47) + (-r48) + (r48_1),
       + (-r50) + (r50_1) + (r53) + (-r53_1),
       + (r53) + (-r53_1) + (-r55) + (r55_1),
       + (r50) + (-r50_1) + (-r52),
       + (r52) + (-r53) + (r53_1),
       + (-r55) + (r55_1) + (r58) + (-r58_1),
       + (r58) + (-r58_1) + (-r60) + (r60_1),
       + (r55) + (-r55_1) + (-r57),
       + (r57) + (-r58) + (r58_1),
       + (-r60) + (r60_1) + (r63) + (-r63_1),
       + (r63) + (-r63_1) + (-r65) + (r65_1),
       + (r60) + (-r60_1) + (-r62),
       + (r62) + (-r63) + (r63_1),
       + (-r65) + (r65_1) + (r68) + (-r68_1),
       + (r68) + (-r68_1) + (-r70) + (r70_1),
       + (r65) + (-r65_1) + (-r67),
       + (r67) + (-r68) + (r68_1),
       + (-r70) + (r70_1) + (r73) + (-r73_1),
       + (r73) + (-r73_1),
       + (r70) + (-r70_1) + (-r72),
       + (r72) + (-r73) + (r73_1),
       + (-r75) + (r75_1) + (r78) + (-r78_1),
       + (r75) + (-r75_1) + (-r77),
       + (r77) + (-r78) + (r78_1)    ])

for i,_ in enumerate(f):
    f[i] = func_array[i]
#func_array[0]

# setup output matrix, only isobutanol x[62]
nouts = 1 
C = np.zeros((nouts, len(x)), dtype = int)
C[0][62] = 1 
C = C.tolist()

# Setup system
x_init_try = np.ones(n)

sys = System(x, f, params = P, params_values = params_values, C = C, x_init = x_init)

sys_ode = get_ODE(sys, timepoints_ode)
sol = sys_ode.solve_system()