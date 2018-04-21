
##############################################
###  Plot MFs in overdensity bins by summing 
###  across all Voronoi redshift slices.
##############################################

plot_fields0 = []
plot_fields0.append(mf_data('N200',    0.691,  0.027))
plot_fields0.append(mf_data('SC1324',  0.755,  0.033))
plot_fields0.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields0.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields0.append(mf_data('N5281',   0.818,  0.029))
plot_fields0.append(mf_data('SC1604',  0.910,  0.029))
plot_fields0.append(mf_data('SC0910',  1.110,  0.035))
plot_fields0.append(mf_data('SC0849',  1.261,  0.029))

cmap = pylab.cm.cool
voronoi_bins = ['_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']
voronoi_labels_merged = ['0.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

voronoi_bins = ['_n05v00', '_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['-0.5 < log(1+$\delta_{\mathrm{gal}}$) < 0.0', '0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

voronoi_bins = ['_n10v05', '_n05v00', '_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['-1.0 < log(1+$\delta_{\mathrm{gal}}$) < -0.5', '-0.5 < log(1+$\delta_{\mathrm{gal}}$) < 0.0', '0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

###  Field Mass Functions - ZFOURGE 2014
###  /Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table1_TOT.dat
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qu = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

###  GCLASS
vdB13 = mypy.readcat('../data/table3_vanderBurg+2013.dat')





if True:

	pylab.ioff()

	for iteration in range(120, 2**len(plot_fields0)):

		iter_name = 'running %3i / %3i :  ' % (iteration, 2**len(plot_fields0))

		###  converting iteration number to binary to indicate which fields to include/exclude
		b = bin(iteration)[2:]
		b = '0'*(8-len(b)) + b

		plot_fields = []
		for i_unit in range(8):
			if b[i_unit] == '1':
				plot_fields.append(plot_fields0[i_unit])
				iter_name += '  %s' % plot_fields0[i_unit].name


		print iter_name





		lmassbins = plot_fields[0].mf_voronoi_slices.lmassbins
		lmassbars = plot_fields[0].mf_voronoi_slices.lmassbars
		dm = lmassbins[1] - lmassbins[0]

		overdens_bins = plot_fields[0].mf_voronoi_slices.overdens_bins
		overdens_bars = plot_fields[0].mf_voronoi_slices.overdens_bars
		d_overdens = overdens_bins[1] - overdens_bins[0]


		data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
		data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
		data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



		###  master 3d array of number counts. Axes are LSSfield, overdensity bin, massbin, 
		master_number_counts = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_number_counts_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_number_counts_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_number_counts_field = pylab.zeros((len(plot_fields), len(lmassbars)))

		master_volumes = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_volumes_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_volumes_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
		master_volumes_field = pylab.zeros((len(plot_fields), len(lmassbars)))

		master_nslices = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))   # total number of voronoi slices

		for fi in range(len(plot_fields)):
			f =  plot_fields[fi]

			###  have to skip by twos because of overlappipng voronoi zslices
			for zi in range(0, len(f.mf_voronoi_slices.zbars), 2):

				###  skipping if zslice is behind cluster
				#if f.zclust < f.mf_voronoi_slices.zlos[zi] - 0.2: continue

				###  skipping redshifts of choice
				#if f.mf_voronoi_slices.zlos[zi] < 1: continue

				lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
				lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

				###  iterating over massbins, checking for mass completeness
				for mi in range(len(lmassbars)):

					if lmassbins[mi] > lmasslimit:
						master_number_counts_field[fi][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
						master_volumes_field[fi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

						master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05[zi][mi]
						master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

						master_number_counts[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00[zi][mi]
						master_volumes[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

						master_number_counts[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05[zi][mi]
						master_volumes[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

						master_number_counts[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10[zi][mi]
						master_volumes[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

						master_number_counts[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15[zi][mi]
						master_volumes[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

						master_number_counts[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20[zi][mi]
						master_volumes[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

						###  counting number of voronoi slices in mass-overdens bins
						if f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi] > 0: master_nslices[fi][0][mi] += 1
						if f.mf_voronoi_slices.volumes_00v05_Mpc3[zi] > 0: master_nslices[fi][1][mi] += 1
						if f.mf_voronoi_slices.volumes_05v10_Mpc3[zi] > 0: master_nslices[fi][2][mi] += 1
						if f.mf_voronoi_slices.volumes_10v15_Mpc3[zi] > 0: master_nslices[fi][3][mi] += 1
						if f.mf_voronoi_slices.volumes_15v20_Mpc3[zi] > 0: master_nslices[fi][4][mi] += 1

					if lmassbins[mi] > lmasslimit:
						master_number_counts_sf[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05_sf[zi][mi]
						master_volumes_sf[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

						master_number_counts_sf[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00_sf[zi][mi]
						master_volumes_sf[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

						master_number_counts_sf[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05_sf[zi][mi]
						master_volumes_sf[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

						master_number_counts_sf[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10_sf[zi][mi]
						master_volumes_sf[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

						master_number_counts_sf[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15_sf[zi][mi]
						master_volumes_sf[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

						master_number_counts_sf[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20_sf[zi][mi]
						master_volumes_sf[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

					if lmassbins[mi] > lmasslimit_ssp:
						master_number_counts_qu[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05_qu[zi][mi]
						master_volumes_qu[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

						master_number_counts_qu[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00_qu[zi][mi]
						master_volumes_qu[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

						master_number_counts_qu[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05_qu[zi][mi]
						master_volumes_qu[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

						master_number_counts_qu[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10_qu[zi][mi]
						master_volumes_qu[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

						master_number_counts_qu[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15_qu[zi][mi]
						master_volumes_qu[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

						master_number_counts_qu[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20_qu[zi][mi]
						master_volumes_qu[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]




		###  plotting results

		fig = pylab.figure(figsize=(15.3, 6.8))
		sp2 = fig.add_subplot(122)
		sp1 = fig.add_subplot(121)

		sp1.grid()
		sp2.grid()
		sp1.minorticks_on()
		sp2.minorticks_on()
		sp1.set_yscale('log')
		sp2.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp2.set_ylabel('Ratio')

		sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])
		sp2.axis([8.4, 11.9, 0.5, 1000])

		fig.subplots_adjust(left=0.09)

		ydatas = []
		fits_tot, covs_tot = [], []
		line_fits_tot, line_covs_tot = [], []
		mlim_hi = 11.55



		###  plotting ZFOURGE
		lmassbars_zf = zfourge_massfunc_tot[:,0]

		area_zf = 316.   # arcmin**2
		volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
		volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
		volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
		volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

		lmf_050z075 = zfourge_massfunc_tot[:,4]
		lmf_075z100 = zfourge_massfunc_tot[:,7]
		lmf_100z125 = zfourge_massfunc_tot[:,10]

		lmf_err_050z075 = (zfourge_massfunc_tot[:,5] + zfourge_massfunc_tot[:,6]) / 2.
		lmf_err_075z100 = (zfourge_massfunc_tot[:,8] + zfourge_massfunc_tot[:,9]) / 2.
		lmf_err_100z125 = (zfourge_massfunc_tot[:,11] + zfourge_massfunc_tot[:,12]) / 2.

		lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
		lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
		lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

		inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))


		sp1.fill_between(lmassbars_zf[inds_zf], 10**(lmf_050z125[inds_zf]-lmf_err_050z125[inds_zf]), 10**(lmf_050z125[inds_zf]+lmf_err_050z125[inds_zf]), color='#cccccc', label='ZFOURGE1')



		###  deciding which overdensity bins to plot
		o_range = range(1, len(overdens_bars))
		offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))

		for vi_counter, vi in enumerate(o_range):
		#for vi in [0]:

			color_i = cmap(vi_counter/(len(o_range)-1.))

			inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
			inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
			x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
			ydatas.append(y)
			ylo, yhi = mypy.massfunc.CI(n)
			elo = y * (ylo / n)
			ehi = y * (yhi / n)


			###  adding to data_table
			data_table_tot[:,0] = lmassbars
			data_table_tot[inds,3*vi_counter+1] = y[inds]
			data_table_tot[inds,3*vi_counter+2] = elo[inds]
			data_table_tot[inds,3*vi_counter+3] = ehi[inds]





		###  writing data_table to ascii file
		outname = '../data/massFunctions_TOT_Tomczak+2017_%ifields' % len(plot_fields)
		for f in plot_fields: outname += '_' + f.name
		outer = open(outname + '.dat', 'w')

		outer.write('#    Galaxy stellar mass functions for all galaxies\n')
		outer.write('#\n')
		outer.write('#    Column  Units            Description\n')
		outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
		outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
		outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
		outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
		outer.write('#\n')
		outer.write('#    Overdensity bins\n')
		outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
		outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
		outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
		outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
		outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
		outer.write('#\n')
		outer.write('# lmass')
		for vi_counter, vi in enumerate(o_range):
			if len(voronoi_bins[vi]) == 6:
				outer.write('   phi%s' % (voronoi_bins[vi]))
				outer.write('   elo%s' % (voronoi_bins[vi]))
				outer.write('   ehi%s' % (voronoi_bins[vi]))
			if len(voronoi_bins[vi]) == 7:
				outer.write('  phi%s' % (voronoi_bins[vi]))
				outer.write('  elo%s' % (voronoi_bins[vi]))
				outer.write('  ehi%s' % (voronoi_bins[vi]))
		outer.write('\n')

		for row in data_table_tot:
			outer.write('%6.2f ' % row[0])
			for item in row[1:]:
				outer.write('  %.4e' % item)
			outer.write('\n')

		outer.close()











		##############################
		###  plotting for SF galaxies
		##############################

		fig = pylab.figure(figsize=(15.5, 12.7))
		sp2 = fig.add_subplot(222)
		sp1 = fig.add_subplot(221)
		sp4 = fig.add_subplot(224)
		sp3 = fig.add_subplot(223)

		#sp1.grid()
		#sp2.grid()
		#sp3.grid()
		#sp4.grid()
		sp1.minorticks_on()
		sp2.minorticks_on()
		sp3.minorticks_on()
		sp4.minorticks_on()
		sp1.set_yscale('log')
		sp2.set_yscale('log')
		sp3.set_yscale('log')
		sp4.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp2.set_ylabel('Ratio')
		sp3.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp4.set_ylabel('Ratio')

		sp1.axis([8.4, 11.9, 8*10**-6, 4*10**-1])
		sp2.axis([8.4, 11.9, 0.5, 2000])
		sp3.axis([8.4, 11.9, 8*10**-6, 4*10**-1])
		sp4.axis([8.4, 11.9, 0.5, 2000])

		fig.subplots_adjust(left=0.09, hspace=0, bottom=0.08)




		###  plotting ZFOURGE
		lmassbars_zf = zfourge_massfunc_sf[:,0]

		area_zf = 316.   # arcmin**2
		volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
		volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
		volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
		volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

		lmf_050z075 = zfourge_massfunc_sf[:,4]
		lmf_075z100 = zfourge_massfunc_sf[:,7]
		lmf_100z125 = zfourge_massfunc_sf[:,10]

		lmf_err_050z075 = (zfourge_massfunc_sf[:,5] + zfourge_massfunc_sf[:,6]) / 2.
		lmf_err_075z100 = (zfourge_massfunc_sf[:,8] + zfourge_massfunc_sf[:,9]) / 2.
		lmf_err_100z125 = (zfourge_massfunc_sf[:,11] + zfourge_massfunc_sf[:,12]) / 2.

		lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
		lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
		lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

		inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))




		ydatas = []
		fits_sf, covs_sf = [], []
		for vi_counter, vi in enumerate(o_range):
		#for vi in [0]:

			color_i = cmap(vi_counter/(len(o_range)-1.))

			inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
			inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
			x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
			ydatas.append(y)
			ylo, yhi = mypy.massfunc.CI(n)

			n_error = n * 1.
			y_error = y * 1.
			ninds = pylab.find(n_error == 0)
			n_error[ninds] = 1.
			y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

			elo = y_error * (ylo / n_error)
			ehi = y_error * (yhi / n_error)

			###  adding to data_table
			data_table_sf[:,0] = lmassbars
			data_table_sf[inds,3*vi_counter+1] = y[inds]
			data_table_sf[inds,3*vi_counter+2] = elo[inds]
			data_table_sf[inds,3*vi_counter+3] = ehi[inds]







		###  writing data_table to ascii file
		outname = '../data/massFunctions_SF_Tomczak+2017_%ifields' % len(plot_fields)
		for f in plot_fields: outname += '_' + f.name
		outer = open(outname + '.dat', 'w')

		outer.write('#    Galaxy stellar mass functions for star-forming galaxies\n')
		outer.write('#\n')
		outer.write('#    Column  Units            Description\n')
		outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
		outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
		outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
		outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
		outer.write('#\n')
		outer.write('#    Overdensity bins\n')
		outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
		outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
		outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
		outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
		outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
		outer.write('#\n')
		outer.write('# lmass')
		for vi_counter, vi in enumerate(o_range):
			if len(voronoi_bins[vi]) == 6:
				outer.write('   phi%s' % (voronoi_bins[vi]))
				outer.write('   elo%s' % (voronoi_bins[vi]))
				outer.write('   ehi%s' % (voronoi_bins[vi]))
			if len(voronoi_bins[vi]) == 7:
				outer.write('  phi%s' % (voronoi_bins[vi]))
				outer.write('  elo%s' % (voronoi_bins[vi]))
				outer.write('  ehi%s' % (voronoi_bins[vi]))
		outer.write('\n')

		for row in data_table_sf:
			outer.write('%6.2f ' % row[0])
			for item in row[1:]:
				outer.write('  %.4e' % item)
			outer.write('\n')

		outer.close()















		##############################
		###  plotting for QU galaxies
		##############################




		###  plotting ZFOURGE
		lmassbars_zf = zfourge_massfunc_qu[:,0]

		area_zf = 316.   # arcmin**2
		volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
		volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
		volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
		volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

		lmf_050z075 = zfourge_massfunc_qu[:,4]
		lmf_075z100 = zfourge_massfunc_qu[:,7]
		lmf_100z125 = zfourge_massfunc_qu[:,10]

		lmf_err_050z075 = (zfourge_massfunc_qu[:,5] + zfourge_massfunc_qu[:,6]) / 2.
		lmf_err_075z100 = (zfourge_massfunc_qu[:,8] + zfourge_massfunc_qu[:,9]) / 2.
		lmf_err_100z125 = (zfourge_massfunc_qu[:,11] + zfourge_massfunc_qu[:,12]) / 2.

		lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
		lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
		lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

		inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))




		qu_inds = []   # indices where qu SMFs are mass-complete for each of the overdensity bins

		ydatas = []
		fits_qu, covs_qu = [], []
		for vi_counter, vi in enumerate(o_range):
		#for vi in [0]:

			color_i = cmap(vi_counter/(len(o_range)-1.))

			#inds = pylab.find((master_volumes_qu.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
			inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
				              (master_volumes_qu.sum(axis=0)[vi] > 0) & 
				              (lmassbars < mlim_hi))
			qu_inds.append(inds)
			x, y, n = lmassbars, (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0) / dm)[vi], master_number_counts_qu.sum(axis=0)[vi]
			ydatas.append(y)
			ylo, yhi = mypy.massfunc.CI(n)
			elo = y * (ylo / n)
			ehi = y * (yhi / n)


			###  adding to data_table
			data_table_qu[:,0] = lmassbars
			data_table_qu[inds,3*vi_counter+1] = y[inds]
			data_table_qu[inds,3*vi_counter+2] = elo[inds]
			data_table_qu[inds,3*vi_counter+3] = ehi[inds]












		###  writing data_table to ascii file
		outname = '../data/massFunctions_QU_Tomczak+2017_%ifields' % len(plot_fields)
		for f in plot_fields: outname += '_' + f.name
		outer = open(outname + '.dat', 'w')

		outer.write('#    Galaxy stellar mass functions for quiescent galaxies\n')
		outer.write('#\n')
		outer.write('#    Column  Units            Description\n')
		outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
		outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
		outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
		outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
		outer.write('#\n')
		outer.write('#    Overdensity bins\n')
		outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
		outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
		outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
		outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
		outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
		outer.write('#\n')
		outer.write('# lmass')
		for vi_counter, vi in enumerate(o_range):
			if len(voronoi_bins[vi]) == 6:
				outer.write('   phi%s' % (voronoi_bins[vi]))
				outer.write('   elo%s' % (voronoi_bins[vi]))
				outer.write('   ehi%s' % (voronoi_bins[vi]))
			if len(voronoi_bins[vi]) == 7:
				outer.write('  phi%s' % (voronoi_bins[vi]))
				outer.write('  elo%s' % (voronoi_bins[vi]))
				outer.write('  ehi%s' % (voronoi_bins[vi]))
		outer.write('\n')

		for row in data_table_qu:
			outer.write('%6.2f ' % row[0])
			for item in row[1:]:
				outer.write('  %.4e' % item)
			outer.write('\n')

		outer.close()

		pylab.close()
		pylab.close()




























############################################################################
###  Generating plots and stats for jack-knife test of the SMF measurements
############################################################################


import glob
import mypy
import pylab
import numpy

def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
	factor1 = pylab.log(10) * pylab.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
	factor2 = phistar1 * 10**(a1*(lmassax - lmstar)) + phistar2 * 10**(a2*(lmassax - lmstar))
	return factor1 * factor2

def schechter_mf(xaxis, alpha, mstar, phistar):
    """
    #  DESCRIPTION:
    #    Returns the values for a Schechter mass function
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis = input mass value(s)
    #    alpha = Schechter parameters
    #    mstar = Schechter parameters
    #    phistar = Schechter parameters    
    """
    return numpy.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * numpy.exp(-10**(xaxis-mstar))


fig = pylab.figure(figsize=(10.5, 9.3))
sp1 = fig.add_subplot(111)

sp1.grid()
sp1.minorticks_on()
sp1.set_yscale('log')

sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')

sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])

fig.subplots_adjust(left=0.1)

mlim_hi = 11.55



files = glob.glob('../data/massFunctions_TOT_Tomczak+2017_[87654321]fields*dat')



schechter_fits = [[] for i_vbin in range(5)]
schechter_covs = [[] for i_vbin in range(5)]
dschechter_fits = [[] for i_vbin in range(5)]
dschechter_covs = [[] for i_vbin in range(5)]
p0_double = [10.9, -0.4,  -1.5,  1.5e-2, 1.2e-3]
p0_single = [-1., 10.9, 5.e-3]


voronoi_labels = ['-0.5 < log(1+$\delta_{\mathrm{gal}}$) < 0.0', '0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']


###  multi-dimensional list to encapsulate all the SMFs from each jack-knife iteration
smfdat_8fields = mypy.readcat('../data/massFunctions_TOT_Tomczak+2017_8fields_N200_SC1324_RCS0224_RXJ1716_N5281_SC1604_SC0910_SC0849.dat')
master_phis = [[[] for i_lmassbin in range(len(smfdat_8fields.lmass))] for i_vbin in range(5)]


widgets = ['  Running: ', progressbar.Percentage(), 
           ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
           ' ', progressbar.ETA(), 
           ' ', progressbar.FileTransferSpeed()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(files)+1).start()

for i_file in range(len(files)):

	pbar.update(i_file)
	data = mypy.readcat(files[i_file])


	phi_n05v00 = data.phi_n05v00
	phi_00v05 = data.phi_00v05
	phi_05v10 = data.phi_05v10
	phi_10v15 = data.phi_10v15
	phi_15v20 = data.phi_15v20

	errs_n05v00 = (data.ehi_n05v00 + data.elo_n05v00) / 2.
	errs_00v05 = (data.ehi_00v05 + data.elo_00v05) / 2.
	errs_05v10 = (data.ehi_05v10 + data.elo_05v10) / 2.
	errs_10v15 = (data.ehi_10v15 + data.elo_10v15) / 2.
	errs_15v20 = (data.ehi_15v20 + data.elo_15v20) / 2.


	phis = [phi_n05v00, phi_00v05, phi_05v10, phi_10v15, phi_15v20]
	errs = [errs_n05v00, errs_00v05, errs_05v10, errs_10v15, errs_15v20]



	for i_vbin in range(len(phis)):

		color_i = pylab.cm.cool(i_vbin / (len(phis) - 1.))

		inds_plot = pylab.find((phis[i_vbin] > 0) & (data.lmass < mlim_hi))

		sp1.plot(data.lmass[inds_plot], phis[i_vbin][inds_plot], color=color_i, lw=2, alpha=0.1, zorder=1)
		if '8fields' in files[i_file]:
			sp1.plot(data.lmass[inds_plot], phis[i_vbin][inds_plot], color='k', lw=1.5, zorder=2)


		###  adding to master phi array
		for i_lmassbin in inds_plot:
			master_phis[i_vbin][i_lmassbin].append(phis[i_vbin][i_lmassbin])




		###  double Schechter fit
		try:
			fit, cov = optimize.curve_fit(dschechter, data.lmass[inds_plot], phis[i_vbin][inds_plot], p0=p0_double, sigma=errs[i_vbin][inds_plot])
			dschechter_fits[i_vbin].append(fit)
			dschechter_covs[i_vbin].append(cov.diagonal()**0.5)

		except:
			pass


		###  single Schechter fit
		try:
			fit, cov = optimize.curve_fit(schechter_mf, data.lmass[inds_plot], phis[i_vbin][inds_plot], p0=p0_single, sigma=errs[i_vbin][inds_plot])
			schechter_fits[i_vbin].append(fit)
			schechter_covs[i_vbin].append(cov.diagonal()**0.5)

		except:
			pass


for i in range(5):
	dschechter_fits[i] = numpy.array(dschechter_fits[i])
	dschechter_covs[i] = numpy.array(dschechter_covs[i])
	schechter_fits[i] = numpy.array(schechter_fits[i])
	schechter_covs[i] = numpy.array(schechter_covs[i])

print '\n\n'




###  print scatter in Schechter parameters
for i in range(5):



	print '%s' % voronoi_labels[i]

	###  single Schechter alpha
	p16 = numpy.percentile(schechter_fits[i][:,0], 16)
	p50 = numpy.percentile(schechter_fits[i][:,0], 50)
	p84 = numpy.percentile(schechter_fits[i][:,0], 84)

	print 'alpha = %5.2f +/- %.2f (%.2f)' % (p50, (p84-p16)/2., mypy.nmad(schechter_fits[i][:,0]))


	###  single Schechter Mstar
	p16 = numpy.percentile(schechter_fits[i][:,1], 16)
	p50 = numpy.percentile(schechter_fits[i][:,1], 50)
	p84 = numpy.percentile(schechter_fits[i][:,1], 84)

	print 'Mstar = %5.2f +/- %.2f (%.2f)' % (p50, (p84-p16)/2., mypy.nmad(schechter_fits[i][:,1]))

	print ''



	###  double Schechter Mstar
	p16 = numpy.percentile(dschechter_fits[i][:,0], 16)
	p50 = numpy.percentile(dschechter_fits[i][:,0], 50)
	p84 = numpy.percentile(dschechter_fits[i][:,0], 84)

	print 'Mstar  = %5.2f +/- %.2f (%.2f)' % (p50, (p84-p16)/2., mypy.nmad(dschechter_fits[i][:,0]))


	###  double Schechter alpha1
	p16 = numpy.percentile(dschechter_fits[i][:,1], 16)
	p50 = numpy.percentile(dschechter_fits[i][:,1], 50)
	p84 = numpy.percentile(dschechter_fits[i][:,1], 84)

	print 'alpha1 = %5.2f +/- %.2f (%.2f)' % (p50, (p84-p16)/2., mypy.nmad(dschechter_fits[i][:,1]))


	###  double Schechter alpha1
	p16 = numpy.percentile(dschechter_fits[i][:,2], 16)
	p50 = numpy.percentile(dschechter_fits[i][:,2], 50)
	p84 = numpy.percentile(dschechter_fits[i][:,2], 84)

	print 'alpha2 = %5.2f +/- %.2f (%.2f)' % (p50, (p84-p16)/2., mypy.nmad(dschechter_fits[i][:,2]))


	print ''
	print ''
	print ''
	print ''
	print ''


















median_matrix  = numpy.zeros((5, len(smfdat_8fields.lmass)))
scatter_matrix = numpy.zeros((5, len(smfdat_8fields.lmass)))

for i_vbin in range(5):
	for i_lmassbin in range(len(smfdat_8fields.lmass)):

		if len(master_phis[i_vbin][i_lmassbin]) > 0:

			p16 = numpy.percentile(master_phis[i_vbin][i_lmassbin], 16)
			p50 = numpy.percentile(master_phis[i_vbin][i_lmassbin], 50)
			p84 = numpy.percentile(master_phis[i_vbin][i_lmassbin], 84)
			
			median_matrix[i_vbin][i_lmassbin] = numpy.log10(p50)
			scatter_matrix[i_vbin][i_lmassbin] = ((p84-p16)/2. / p50) / numpy.log(10)

			sp1.plot(smfdat_8fields.lmass[i_lmassbin], p50, 'ko')



inds = smfdat_8fields.phi_n05v00 > 0
print 'max offset at -0.5<log(1+d)<0.0 =', abs(numpy.log10(smfdat_8fields.phi_n05v00[inds]) - median_matrix[0][inds])

inds = smfdat_8fields.phi_00v05 > 0
print 'max offset at  0.0<log(1+d)<0.5 =', abs(numpy.log10(smfdat_8fields.phi_00v05[inds]) - median_matrix[1][inds]).max()

inds = smfdat_8fields.phi_05v10 > 0
print 'max offset at  0.5<log(1+d)<1.0 =', abs(numpy.log10(smfdat_8fields.phi_05v10[inds]) - median_matrix[2][inds]).max()

inds = smfdat_8fields.phi_10v15 > 0
print 'max offset at  1.0<log(1+d)<1.5 =', abs(numpy.log10(smfdat_8fields.phi_10v15[inds]) - median_matrix[3][inds]).max()

inds = smfdat_8fields.phi_15v20 > 0
print 'max offset at  1.5<log(1+d)<2.0 =', abs(numpy.log10(smfdat_8fields.phi_15v20[inds]) - median_matrix[4][inds]).max()


print ''


s1 = numpy.sort(scatter_matrix[0][scatter_matrix[0] > 0])
s2 = numpy.sort(scatter_matrix[1][scatter_matrix[1] > 0])
s3 = numpy.sort(scatter_matrix[2][scatter_matrix[2] > 0])
s4 = numpy.sort(scatter_matrix[3][scatter_matrix[3] > 0])
s5 = numpy.sort(scatter_matrix[4][scatter_matrix[4] > 0])

print 'max scatter at -0.5<log(1+d)<0.0 =', s1[-4:]
print 'max scatter at  0.0<log(1+d)<0.5 =', s2[-4:]
print 'max scatter at  0.5<log(1+d)<1.0 =', s3[-4:]
print 'max scatter at  1.0<log(1+d)<1.5 =', s4[-4:]
print 'max scatter at  1.5<log(1+d)<2.0 =', s5[-4:]

print ''

print 'median scatter at -0.5<log(1+d)<0.0 =', numpy.median(s1)
print 'median scatter at  0.0<log(1+d)<0.5 =', numpy.median(s2)
print 'median scatter at  0.5<log(1+d)<1.0 =', numpy.median(s3)
print 'median scatter at  1.0<log(1+d)<1.5 =', numpy.median(s4)
print 'median scatter at  1.5<log(1+d)<2.0 =', numpy.median(s5)









