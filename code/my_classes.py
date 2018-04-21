
import mypy
import numpy




class merger_simulation:
	def __init__(self, lmassax):
		self.lmassax = lmassax
		self.dm = abs(lmassax[1] - lmassax[0])
		self.fraction_merged = numpy.array([])
		self.ngals_fmerged = []                # SMFs at each value in fraction_merged

		self.n_initial = 0
		self.mass_initial = 0.

		self.n_initial_gt_9 = 0
		self.n_merged_gt_9 = 0

		self.mass_in_icm = numpy.array([])
		self.f_icm = numpy.array([])
		self.ngals_icm = []                 # SMFs at each value in mass_in_icl

		self.n_vminor_mergers = []           # mean number of very minor mergers in bins of lmass (mu < 1:10)
		self.n_minor_mergers = []           # mean number of minor mergers in bins of lmass (1:10 < mu < 1:4)
		self.n_major_mergers = []           # mean number of major mergers in bins of lmass (mu > 1:4)






