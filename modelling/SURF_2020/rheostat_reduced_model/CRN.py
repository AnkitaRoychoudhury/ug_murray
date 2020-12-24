import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class SBMLmodel(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters

		self.c = {} #Dictionary of compartments
		self.c['default'] = sbmltoodepy.modelclasses.Compartment(1e-06, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("default"))

		self.s = {} #Dictionary of chemical species
		self.s['enzyme_hex'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_hex"))
		self.s['metabolite_atp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("metabolite_atp"))
		self.s['molecule_glucose'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_glucose"))
		self.s['molecule_g6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_g6p"))
		self.s['metabolite_adp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("metabolite_adp"))
		self.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_hex_metabolite_atp_molecule_glucose"))
		self.s['complex_enzyme_hex_metabolite_adp_molecule_g6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_hex_metabolite_adp_molecule_g6p"))
		self.s['enzyme_pgi'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_pgi"))
		self.s['molecule_f6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_f6p"))
		self.s['complex_enzyme_pgi_molecule_g6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgi_molecule_g6p"))
		self.s['complex_enzyme_pgi_molecule_f6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgi_molecule_f6p"))
		self.s['enzyme_pfk'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_pfk"))
		self.s['molecule_f16p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_f16p"))
		self.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pfk_metabolite_atp_molecule_f6p"))
		self.s['complex_enzyme_pfk_metabolite_adp_molecule_f16p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pfk_metabolite_adp_molecule_f16p"))
		self.s['enzyme_ald_tpi'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_ald_tpi"))
		self.s['molecule_g3p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_g3p"))
		self.s['complex_enzyme_ald_tpi_molecule_f16p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_ald_tpi_molecule_f16p"))
		self.s['complex_enzyme_ald_tpi_2x_molecule_g3p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_ald_tpi_2x_molecule_g3p"))
		self.s['enzyme_gapN'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_gapN"))
		self.s['metabolite_nadp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("metabolite_nadp"))
		self.s['molecule_3pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_3pg"))
		self.s['metabolite_nadph'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("metabolite_nadph"))
		self.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p"))
		self.s['complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg"))
		self.s['enzyme_mGapDH'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_mGapDH"))
		self.s['metabolite_pi'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("metabolite_pi"))
		self.s['molecule_13bpg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_13bpg"))
		self.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p"))
		self.s['complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg"))
		self.s['enzyme_pgk'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_pgk"))
		self.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgk_metabolite_adp_molecule_13bpg"))
		self.s['complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg"))
		self.s['enzyme_pgm'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_pgm"))
		self.s['molecule_2pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_2pg"))
		self.s['complex_enzyme_pgm_2x_molecule_3pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgm_2x_molecule_3pg"))
		self.s['complex_enzyme_pgm_2x_molecule_2pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pgm_2x_molecule_2pg"))
		self.s['enzyme_eno'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_eno"))
		self.s['molecule_pep'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_pep"))
		self.s['complex_enzyme_eno_2x_molecule_2pg'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_eno_2x_molecule_2pg"))
		self.s['complex_enzyme_eno_2x_molecule_pep'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_eno_2x_molecule_pep"))
		self.s['enzyme_pyk'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_pyk"))
		self.s['molecule_pyruvate'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_pyruvate"))
		self.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep"))
		self.s['complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate"))
		self.s['enzyme_alsS'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_alsS"))
		self.s['molecule_acetolac'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_acetolac"))
		self.s['complex_enzyme_alsS_2x_molecule_pyruvate'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_alsS_2x_molecule_pyruvate"))
		self.s['complex_enzyme_alsS_molecule_acetolac'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_alsS_molecule_acetolac"))
		self.s['enzyme_IlvC'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_IlvC"))
		self.s['molecule_23dih3mebut'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_23dih3mebut"))
		self.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac"))
		self.s['complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut"))
		self.s['enzyme_IlvD'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_IlvD"))
		self.s['molecule_3me2oxo'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_3me2oxo"))
		self.s['complex_enzyme_IlvD_molecule_23dih3mebut'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_IlvD_molecule_23dih3mebut"))
		self.s['complex_enzyme_IlvD_molecule_3me2oxo'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_IlvD_molecule_3me2oxo"))
		self.s['enzyme_kivD'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_kivD"))
		self.s['molecule_isobutanal'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_isobutanal"))
		self.s['complex_enzyme_kivD_molecule_3me2oxo'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_kivD_molecule_3me2oxo"))
		self.s['complex_enzyme_kivD_molecule_isobutanal'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_kivD_molecule_isobutanal"))
		self.s['enzyme_yahk'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_yahk"))
		self.s['molecule_isobutanol'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("molecule_isobutanol"))
		self.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal"))
		self.s['complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol"))
		self.s['enzyme_atpase'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("enzyme_atpase"))
		self.s['complex_enzyme_atpase_metabolite_atp'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_atpase_metabolite_atp"))
		self.s['complex_enzyme_atpase_metabolite_adp_metabolite_pi'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['default'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("complex_enzyme_atpase_metabolite_adp_metabolite_pi"))

		self.r = {} #Dictionary of reactions
		self.r['r0'] = r0(self)
		self.r['r0_1'] = r0_1(self)
		self.r['r2'] = r2(self)
		self.r['r3'] = r3(self)
		self.r['r3_1'] = r3_1(self)
		self.r['r5'] = r5(self)
		self.r['r5_1'] = r5_1(self)
		self.r['r7'] = r7(self)
		self.r['r8'] = r8(self)
		self.r['r8_1'] = r8_1(self)
		self.r['r10'] = r10(self)
		self.r['r10_1'] = r10_1(self)
		self.r['r12'] = r12(self)
		self.r['r13'] = r13(self)
		self.r['r13_1'] = r13_1(self)
		self.r['r15'] = r15(self)
		self.r['r15_1'] = r15_1(self)
		self.r['r17'] = r17(self)
		self.r['r18'] = r18(self)
		self.r['r18_1'] = r18_1(self)
		self.r['r20'] = r20(self)
		self.r['r20_1'] = r20_1(self)
		self.r['r22'] = r22(self)
		self.r['r23'] = r23(self)
		self.r['r23_1'] = r23_1(self)
		self.r['r25'] = r25(self)
		self.r['r25_1'] = r25_1(self)
		self.r['r27'] = r27(self)
		self.r['r28'] = r28(self)
		self.r['r28_1'] = r28_1(self)
		self.r['r30'] = r30(self)
		self.r['r30_1'] = r30_1(self)
		self.r['r32'] = r32(self)
		self.r['r33'] = r33(self)
		self.r['r33_1'] = r33_1(self)
		self.r['r35'] = r35(self)
		self.r['r35_1'] = r35_1(self)
		self.r['r37'] = r37(self)
		self.r['r38'] = r38(self)
		self.r['r38_1'] = r38_1(self)
		self.r['r40'] = r40(self)
		self.r['r40_1'] = r40_1(self)
		self.r['r42'] = r42(self)
		self.r['r43'] = r43(self)
		self.r['r43_1'] = r43_1(self)
		self.r['r45'] = r45(self)
		self.r['r45_1'] = r45_1(self)
		self.r['r47'] = r47(self)
		self.r['r48'] = r48(self)
		self.r['r48_1'] = r48_1(self)
		self.r['r50'] = r50(self)
		self.r['r50_1'] = r50_1(self)
		self.r['r52'] = r52(self)
		self.r['r53'] = r53(self)
		self.r['r53_1'] = r53_1(self)
		self.r['r55'] = r55(self)
		self.r['r55_1'] = r55_1(self)
		self.r['r57'] = r57(self)
		self.r['r58'] = r58(self)
		self.r['r58_1'] = r58_1(self)
		self.r['r60'] = r60(self)
		self.r['r60_1'] = r60_1(self)
		self.r['r62'] = r62(self)
		self.r['r63'] = r63(self)
		self.r['r63_1'] = r63_1(self)
		self.r['r65'] = r65(self)
		self.r['r65_1'] = r65_1(self)
		self.r['r67'] = r67(self)
		self.r['r68'] = r68(self)
		self.r['r68_1'] = r68_1(self)
		self.r['r70'] = r70(self)
		self.r['r70_1'] = r70_1(self)
		self.r['r72'] = r72(self)
		self.r['r73'] = r73(self)
		self.r['r73_1'] = r73_1(self)
		self.r['r75'] = r75(self)
		self.r['r75_1'] = r75_1(self)
		self.r['r77'] = r77(self)
		self.r['r78'] = r78(self)
		self.r['r78_1'] = r78_1(self)

		self.f = {} #Dictionary of function definitions
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['enzyme_hex'].amount, self.s['metabolite_atp'].amount, self.s['molecule_glucose'].amount, self.s['molecule_g6p'].amount, self.s['metabolite_adp'].amount, self.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'].amount, self.s['complex_enzyme_hex_metabolite_adp_molecule_g6p'].amount, self.s['enzyme_pgi'].amount, self.s['molecule_f6p'].amount, self.s['complex_enzyme_pgi_molecule_g6p'].amount, self.s['complex_enzyme_pgi_molecule_f6p'].amount, self.s['enzyme_pfk'].amount, self.s['molecule_f16p'].amount, self.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'].amount, self.s['complex_enzyme_pfk_metabolite_adp_molecule_f16p'].amount, self.s['enzyme_ald_tpi'].amount, self.s['molecule_g3p'].amount, self.s['complex_enzyme_ald_tpi_molecule_f16p'].amount, self.s['complex_enzyme_ald_tpi_2x_molecule_g3p'].amount, self.s['enzyme_gapN'].amount, self.s['metabolite_nadp'].amount, self.s['molecule_3pg'].amount, self.s['metabolite_nadph'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg'].amount, self.s['enzyme_mGapDH'].amount, self.s['metabolite_pi'].amount, self.s['molecule_13bpg'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg'].amount, self.s['enzyme_pgk'].amount, self.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'].amount, self.s['complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg'].amount, self.s['enzyme_pgm'].amount, self.s['molecule_2pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_3pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_2pg'].amount, self.s['enzyme_eno'].amount, self.s['molecule_pep'].amount, self.s['complex_enzyme_eno_2x_molecule_2pg'].amount, self.s['complex_enzyme_eno_2x_molecule_pep'].amount, self.s['enzyme_pyk'].amount, self.s['molecule_pyruvate'].amount, self.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'].amount, self.s['complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate'].amount, self.s['enzyme_alsS'].amount, self.s['molecule_acetolac'].amount, self.s['complex_enzyme_alsS_2x_molecule_pyruvate'].amount, self.s['complex_enzyme_alsS_molecule_acetolac'].amount, self.s['enzyme_IlvC'].amount, self.s['molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'].amount, self.s['complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut'].amount, self.s['enzyme_IlvD'].amount, self.s['molecule_3me2oxo'].amount, self.s['complex_enzyme_IlvD_molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvD_molecule_3me2oxo'].amount, self.s['enzyme_kivD'].amount, self.s['molecule_isobutanal'].amount, self.s['complex_enzyme_kivD_molecule_3me2oxo'].amount, self.s['complex_enzyme_kivD_molecule_isobutanal'].amount, self.s['enzyme_yahk'].amount, self.s['molecule_isobutanol'].amount, self.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'].amount, self.s['complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol'].amount, self.s['enzyme_atpase'].amount, self.s['complex_enzyme_atpase_metabolite_atp'].amount, self.s['complex_enzyme_atpase_metabolite_adp_metabolite_pi'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[-1,1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[-1,1,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,-1,1,0,0,0.],[-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,1,-1,0,0,0,0,0,0,0,0,1,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,-2,2,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1.],[ 1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,-1,1,0,1,-1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,1,-1,1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,-1,1,0,1,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0.,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,2,-2,-2,2,0,0,0,-2,2,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-1,1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-2,2,0,0,0,-2,2,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,-1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,2,-2,0,0,0,0,0,0,0,0,2,-2,-2.,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,2,-2,0,0,0,2,-2,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,-1,1,0,1,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,-1,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,1,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.,1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,-1,1,0,1,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,2,-2,-2,2,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1,1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,-1,1,0,1,-1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,-2,2,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,1.,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1.,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,-1,1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,-1,1,0,1,-1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,1,-1,1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,1,-1,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,0,0,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,0,0.,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1.,0,1,-1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,-1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1.,-1,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,-1,1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,-1,1,0,1,-1.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,-1,-1,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,1,-1,1.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['r0'](), self.r['r0_1'](), self.r['r2'](), self.r['r3'](), self.r['r3_1'](), self.r['r5'](), self.r['r5_1'](), self.r['r7'](), self.r['r8'](), self.r['r8_1'](), self.r['r10'](), self.r['r10_1'](), self.r['r12'](), self.r['r13'](), self.r['r13_1'](), self.r['r15'](), self.r['r15_1'](), self.r['r17'](), self.r['r18'](), self.r['r18_1'](), self.r['r20'](), self.r['r20_1'](), self.r['r22'](), self.r['r23'](), self.r['r23_1'](), self.r['r25'](), self.r['r25_1'](), self.r['r27'](), self.r['r28'](), self.r['r28_1'](), self.r['r30'](), self.r['r30_1'](), self.r['r32'](), self.r['r33'](), self.r['r33_1'](), self.r['r35'](), self.r['r35_1'](), self.r['r37'](), self.r['r38'](), self.r['r38_1'](), self.r['r40'](), self.r['r40_1'](), self.r['r42'](), self.r['r43'](), self.r['r43_1'](), self.r['r45'](), self.r['r45_1'](), self.r['r47'](), self.r['r48'](), self.r['r48_1'](), self.r['r50'](), self.r['r50_1'](), self.r['r52'](), self.r['r53'](), self.r['r53_1'](), self.r['r55'](), self.r['r55_1'](), self.r['r57'](), self.r['r58'](), self.r['r58_1'](), self.r['r60'](), self.r['r60_1'](), self.r['r62'](), self.r['r63'](), self.r['r63_1'](), self.r['r65'](), self.r['r65_1'](), self.r['r67'](), self.r['r68'](), self.r['r68_1'](), self.r['r70'](), self.r['r70_1'](), self.r['r72'](), self.r['r73'](), self.r['r73_1'](), self.r['r75'](), self.r['r75_1'](), self.r['r77'](), self.r['r78'](), self.r['r78_1']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['enzyme_hex'].amount, self.s['metabolite_atp'].amount, self.s['molecule_glucose'].amount, self.s['molecule_g6p'].amount, self.s['metabolite_adp'].amount, self.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'].amount, self.s['complex_enzyme_hex_metabolite_adp_molecule_g6p'].amount, self.s['enzyme_pgi'].amount, self.s['molecule_f6p'].amount, self.s['complex_enzyme_pgi_molecule_g6p'].amount, self.s['complex_enzyme_pgi_molecule_f6p'].amount, self.s['enzyme_pfk'].amount, self.s['molecule_f16p'].amount, self.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'].amount, self.s['complex_enzyme_pfk_metabolite_adp_molecule_f16p'].amount, self.s['enzyme_ald_tpi'].amount, self.s['molecule_g3p'].amount, self.s['complex_enzyme_ald_tpi_molecule_f16p'].amount, self.s['complex_enzyme_ald_tpi_2x_molecule_g3p'].amount, self.s['enzyme_gapN'].amount, self.s['metabolite_nadp'].amount, self.s['molecule_3pg'].amount, self.s['metabolite_nadph'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg'].amount, self.s['enzyme_mGapDH'].amount, self.s['metabolite_pi'].amount, self.s['molecule_13bpg'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg'].amount, self.s['enzyme_pgk'].amount, self.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'].amount, self.s['complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg'].amount, self.s['enzyme_pgm'].amount, self.s['molecule_2pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_3pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_2pg'].amount, self.s['enzyme_eno'].amount, self.s['molecule_pep'].amount, self.s['complex_enzyme_eno_2x_molecule_2pg'].amount, self.s['complex_enzyme_eno_2x_molecule_pep'].amount, self.s['enzyme_pyk'].amount, self.s['molecule_pyruvate'].amount, self.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'].amount, self.s['complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate'].amount, self.s['enzyme_alsS'].amount, self.s['molecule_acetolac'].amount, self.s['complex_enzyme_alsS_2x_molecule_pyruvate'].amount, self.s['complex_enzyme_alsS_molecule_acetolac'].amount, self.s['enzyme_IlvC'].amount, self.s['molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'].amount, self.s['complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut'].amount, self.s['enzyme_IlvD'].amount, self.s['molecule_3me2oxo'].amount, self.s['complex_enzyme_IlvD_molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvD_molecule_3me2oxo'].amount, self.s['enzyme_kivD'].amount, self.s['molecule_isobutanal'].amount, self.s['complex_enzyme_kivD_molecule_3me2oxo'].amount, self.s['complex_enzyme_kivD_molecule_isobutanal'].amount, self.s['enzyme_yahk'].amount, self.s['molecule_isobutanol'].amount, self.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'].amount, self.s['complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol'].amount, self.s['enzyme_atpase'].amount, self.s['complex_enzyme_atpase_metabolite_atp'].amount, self.s['complex_enzyme_atpase_metabolite_adp_metabolite_pi'].amount], dtype = np.float64)
		self.s['enzyme_hex'].amount, self.s['metabolite_atp'].amount, self.s['molecule_glucose'].amount, self.s['molecule_g6p'].amount, self.s['metabolite_adp'].amount, self.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'].amount, self.s['complex_enzyme_hex_metabolite_adp_molecule_g6p'].amount, self.s['enzyme_pgi'].amount, self.s['molecule_f6p'].amount, self.s['complex_enzyme_pgi_molecule_g6p'].amount, self.s['complex_enzyme_pgi_molecule_f6p'].amount, self.s['enzyme_pfk'].amount, self.s['molecule_f16p'].amount, self.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'].amount, self.s['complex_enzyme_pfk_metabolite_adp_molecule_f16p'].amount, self.s['enzyme_ald_tpi'].amount, self.s['molecule_g3p'].amount, self.s['complex_enzyme_ald_tpi_molecule_f16p'].amount, self.s['complex_enzyme_ald_tpi_2x_molecule_g3p'].amount, self.s['enzyme_gapN'].amount, self.s['metabolite_nadp'].amount, self.s['molecule_3pg'].amount, self.s['metabolite_nadph'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'].amount, self.s['complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg'].amount, self.s['enzyme_mGapDH'].amount, self.s['metabolite_pi'].amount, self.s['molecule_13bpg'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'].amount, self.s['complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg'].amount, self.s['enzyme_pgk'].amount, self.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'].amount, self.s['complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg'].amount, self.s['enzyme_pgm'].amount, self.s['molecule_2pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_3pg'].amount, self.s['complex_enzyme_pgm_2x_molecule_2pg'].amount, self.s['enzyme_eno'].amount, self.s['molecule_pep'].amount, self.s['complex_enzyme_eno_2x_molecule_2pg'].amount, self.s['complex_enzyme_eno_2x_molecule_pep'].amount, self.s['enzyme_pyk'].amount, self.s['molecule_pyruvate'].amount, self.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'].amount, self.s['complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate'].amount, self.s['enzyme_alsS'].amount, self.s['molecule_acetolac'].amount, self.s['complex_enzyme_alsS_2x_molecule_pyruvate'].amount, self.s['complex_enzyme_alsS_molecule_acetolac'].amount, self.s['enzyme_IlvC'].amount, self.s['molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'].amount, self.s['complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut'].amount, self.s['enzyme_IlvD'].amount, self.s['molecule_3me2oxo'].amount, self.s['complex_enzyme_IlvD_molecule_23dih3mebut'].amount, self.s['complex_enzyme_IlvD_molecule_3me2oxo'].amount, self.s['enzyme_kivD'].amount, self.s['molecule_isobutanal'].amount, self.s['complex_enzyme_kivD_molecule_3me2oxo'].amount, self.s['complex_enzyme_kivD_molecule_isobutanal'].amount, self.s['enzyme_yahk'].amount, self.s['molecule_isobutanol'].amount, self.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'].amount, self.s['complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol'].amount, self.s['enzyme_atpase'].amount, self.s['complex_enzyme_atpase_metabolite_atp'].amount, self.s['complex_enzyme_atpase_metabolite_adp_metabolite_pi'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class r0:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r0")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_hex'].concentration * self.parent.s['metabolite_atp'].concentration * self.parent.s['molecule_glucose'].concentration

class r0_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r0_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'].concentration

class r2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r2")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_hex_metabolite_atp_molecule_glucose'].concentration

class r3:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r3")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_hex_metabolite_adp_molecule_g6p'].concentration

class r3_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r3_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_hex'].concentration * self.parent.s['molecule_g6p'].concentration * self.parent.s['metabolite_adp'].concentration

class r5:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r5")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgi'].concentration * self.parent.s['molecule_g6p'].concentration

class r5_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r5_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgi_molecule_g6p'].concentration

class r7:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r7")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgi_molecule_g6p'].concentration

class r8:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r8")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgi_molecule_f6p'].concentration

class r8_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r8_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgi'].concentration * self.parent.s['molecule_f6p'].concentration

class r10:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r10")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pfk'].concentration * self.parent.s['metabolite_atp'].concentration * self.parent.s['molecule_f6p'].concentration

class r10_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r10_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'].concentration

class r12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r12")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pfk_metabolite_atp_molecule_f6p'].concentration

class r13:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r13")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pfk_metabolite_adp_molecule_f16p'].concentration

class r13_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r13_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pfk'].concentration * self.parent.s['molecule_f16p'].concentration * self.parent.s['metabolite_adp'].concentration

class r15:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r15")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_ald_tpi'].concentration * self.parent.s['molecule_f16p'].concentration

class r15_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r15_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_ald_tpi_molecule_f16p'].concentration

class r17:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r17")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_ald_tpi_molecule_f16p'].concentration

class r18:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r18")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_ald_tpi_2x_molecule_g3p'].concentration

class r18_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r18_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_ald_tpi'].concentration * self.parent.s['molecule_g3p'].concentration**2

class r20:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r20")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_gapN'].concentration * self.parent.s['metabolite_nadp'].concentration**2 * self.parent.s['molecule_g3p'].concentration**2

class r20_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r20_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'].concentration

class r22:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r22")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_gapN_2x_metabolite_nadp_2x_molecule_g3p'].concentration

class r23:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r23")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_gapN_2x_metabolite_nadph_2x_molecule_3pg'].concentration

class r23_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r23_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_gapN'].concentration * self.parent.s['molecule_3pg'].concentration**2 * self.parent.s['metabolite_nadph'].concentration**2

class r25:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r25")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_mGapDH'].concentration * self.parent.s['metabolite_pi'].concentration * self.parent.s['metabolite_nadp'].concentration**2 * self.parent.s['molecule_g3p'].concentration**2

class r25_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r25_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'].concentration

class r27:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r27")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_mGapDH_2x_metabolite_nadp_metabolite_pi_2x_molecule_g3p'].concentration

class r28:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r28")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_mGapDH_2x_metabolite_nadph_molecule_13bpg'].concentration

class r28_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r28_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_mGapDH'].concentration * self.parent.s['molecule_13bpg'].concentration * self.parent.s['metabolite_nadph'].concentration**2

class r30:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r30")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgk'].concentration * self.parent.s['metabolite_adp'].concentration * self.parent.s['molecule_13bpg'].concentration

class r30_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r30_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'].concentration

class r32:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r32")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgk_metabolite_adp_molecule_13bpg'].concentration

class r33:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r33")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgk_metabolite_atp_2x_molecule_3pg'].concentration

class r33_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r33_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgk'].concentration * self.parent.s['molecule_3pg'].concentration**2 * self.parent.s['metabolite_atp'].concentration

class r35:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r35")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgm'].concentration * self.parent.s['molecule_3pg'].concentration**2

class r35_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r35_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgm_2x_molecule_3pg'].concentration

class r37:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r37")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgm_2x_molecule_3pg'].concentration

class r38:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r38")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pgm_2x_molecule_2pg'].concentration

class r38_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r38_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pgm'].concentration * self.parent.s['molecule_2pg'].concentration**2

class r40:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r40")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_eno'].concentration * self.parent.s['molecule_2pg'].concentration**2

class r40_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r40_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_eno_2x_molecule_2pg'].concentration

class r42:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r42")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_eno_2x_molecule_2pg'].concentration

class r43:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r43")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_eno_2x_molecule_pep'].concentration

class r43_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r43_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_eno'].concentration * self.parent.s['molecule_pep'].concentration**2

class r45:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r45")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pyk'].concentration * self.parent.s['metabolite_adp'].concentration**2 * self.parent.s['molecule_pep'].concentration**2

class r45_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r45_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'].concentration

class r47:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r47")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pyk_2x_metabolite_adp_2x_molecule_pep'].concentration

class r48:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r48")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_pyk_2x_metabolite_atp_2x_molecule_pyruvate'].concentration

class r48_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r48_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_pyk'].concentration * self.parent.s['molecule_pyruvate'].concentration**2 * self.parent.s['metabolite_atp'].concentration**2

class r50:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r50")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_alsS'].concentration * self.parent.s['molecule_pyruvate'].concentration**2

class r50_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r50_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_alsS_2x_molecule_pyruvate'].concentration

class r52:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r52")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_alsS_2x_molecule_pyruvate'].concentration

class r53:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r53")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_alsS_molecule_acetolac'].concentration

class r53_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r53_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_alsS'].concentration * self.parent.s['molecule_acetolac'].concentration

class r55:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r55")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_IlvC'].concentration * self.parent.s['metabolite_nadph'].concentration * self.parent.s['molecule_acetolac'].concentration

class r55_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r55_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'].concentration

class r57:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r57")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvC_metabolite_nadph_molecule_acetolac'].concentration

class r58:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r58")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvC_metabolite_nadp_molecule_23dih3mebut'].concentration

class r58_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r58_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_IlvC'].concentration * self.parent.s['molecule_23dih3mebut'].concentration * self.parent.s['metabolite_nadp'].concentration

class r60:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r60")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_IlvD'].concentration * self.parent.s['molecule_23dih3mebut'].concentration

class r60_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r60_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvD_molecule_23dih3mebut'].concentration

class r62:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r62")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvD_molecule_23dih3mebut'].concentration

class r63:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r63")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_IlvD_molecule_3me2oxo'].concentration

class r63_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r63_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_IlvD'].concentration * self.parent.s['molecule_3me2oxo'].concentration

class r65:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r65")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_kivD'].concentration * self.parent.s['molecule_3me2oxo'].concentration

class r65_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r65_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_kivD_molecule_3me2oxo'].concentration

class r67:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r67")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_kivD_molecule_3me2oxo'].concentration

class r68:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r68")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_kivD_molecule_isobutanal'].concentration

class r68_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r68_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_kivD'].concentration * self.parent.s['molecule_isobutanal'].concentration

class r70:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r70")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_yahk'].concentration * self.parent.s['metabolite_nadph'].concentration * self.parent.s['molecule_isobutanal'].concentration

class r70_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r70_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'].concentration

class r72:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r72")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(36000.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_yahk_metabolite_nadph_molecule_isobutanal'].concentration

class r73:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r73")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_yahk_metabolite_nadp_molecule_isobutanol'].concentration

class r73_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r73_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_yahk'].concentration * self.parent.s['molecule_isobutanol'].concentration * self.parent.s['metabolite_nadp'].concentration

class r75:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r75")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_atpase'].concentration * self.parent.s['metabolite_atp'].concentration

class r75_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r75_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_atpase_metabolite_atp'].concentration

class r77:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r77")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(2.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_atpase_metabolite_atp'].concentration

class r78:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r78")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(15.0, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['complex_enzyme_atpase_metabolite_adp_metabolite_pi'].concentration

class r78_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("r78_1")
		self.p['k'] = sbmltoodepy.modelclasses.Parameter(1.5, 'k')

	def __call__(self):
		return self.p['k'].value * self.parent.s['enzyme_atpase'].concentration * self.parent.s['metabolite_adp'].concentration * self.parent.s['metabolite_pi'].concentration

