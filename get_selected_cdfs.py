from pylab import *

def geom(vals):
	val = vals[0]
	for i in range(1,len(vals)):
		val *= vals[i]
	return val**(1./len(vals))

# CMB 
cdf_CMB = loadtxt('../../new/pd_elastic_new_CDF_CMB.txt')
file_max = open('pd_elastic_CDF_CMB_max.txt','w')
file_min = open('pd_elastic_CDF_CMB_min.txt','w')
file_mean = open('pd_elastic_CDF_CMB_mean.txt','w')
file_geom = open('pd_elastic_CDF_CMB_geom.txt','w')
hdr_min = '# cumulative interaction rate CDF with CMB for minimum eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_max = '# cumulative interaction rate CDF with CMB for maximum eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_mean = '# cumulative interaction rate CDF with CMB for mean eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_geom = '# cumulative interaction rate CDF with CMB for geometric mean eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
file_min.write(hdr_min)
file_max.write(hdr_max)
file_mean.write(hdr_mean)
file_geom.write(hdr_geom)

gamma = unique(cdf_CMB[:,2])

for gam in gamma:
	idx = gam==cdf_CMB[:,2]
	dat = cdf_CMB[:,3:][idx]
	max_eps_cdf = zeros(len(dat[0,:]))
	min_eps_cdf = zeros(len(dat[0,:]))
	mean_eps_cdf = zeros(len(dat[0,:]))
	geom_eps_cdf = zeros(len(dat[0,:]))
	for i in range(0,len(dat[0,:])):
		max_eps_cdf[i] = min(dat[:,i])
		min_eps_cdf[i] = max(dat[:,i])
		mean_eps_cdf[i] = mean(dat[:,i])
		geom_eps_cdf[i] = geom(dat[:,i])
	fmt = '%g' + '\t%g'*len(dat[0,:]) + '\n'
	file_max.write(fmt % ((gam,)+tuple(max_eps_cdf)))
	file_min.write(fmt % ((gam,)+tuple(min_eps_cdf)))
	file_mean.write(fmt % ((gam,)+tuple(mean_eps_cdf)))
	file_geom.write(fmt % ((gam,)+tuple(geom_eps_cdf)))
file_max.close()
file_min.close()
file_mean.close()
file_geom.close()

# IRB
cdf_IRB = loadtxt('../../new/pd_elastic_new_CDF_IRB_Kneiske04.txt')
file_max = open('pd_elastic_CDF_IRB_Kneiske04_max.txt','w')
file_min = open('pd_elastic_CDF_IRB_Kneiske04_min.txt','w')
file_mean = open('pd_elastic_CDF_IRB_Kneiske04_mean.txt','w')
file_geom = open('pd_elastic_CDF_IRB_Kneiske04_geom.txt','w')
hdr_min = '# cumulative interaction rate CDF with IRB_Kneiske04 for minimum eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_max = '# cumulative interaction rate CDF with IRB_Kneiske04 for maximum eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_mean = '# cumulative interaction rate CDF with IRB_Kneiske04 for mean eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
hdr_geom = '# cumulative interaction rate CDF with IRB_Kneiske04 for geometric mean eps values\n# log10(gamma)\tcdf for eps from 0.1 MeV to 26.3 GeV in 513 steps\n' 
file_min.write(hdr_min)
file_max.write(hdr_max)
file_mean.write(hdr_mean)
file_geom.write(hdr_geom)

gamma = unique(cdf_IRB[:,2])

for gam in gamma:
	idx = gam==cdf_IRB[:,2]
	dat = cdf_IRB[:,3:][idx]
	max_eps_cdf = zeros(len(dat[0,:]))
	min_eps_cdf = zeros(len(dat[0,:]))
	mean_eps_cdf = zeros(len(dat[0,:]))
	geom_eps_cdf = zeros(len(dat[0,:]))
	for i in range(0,len(dat[0,:])):
		max_eps_cdf[i] = min(dat[:,i])
		min_eps_cdf[i] = max(dat[:,i])
		mean_eps_cdf[i] = mean(dat[:,i])
		geom_eps_cdf[i] = geom(dat[:,i])
	fmt = '%g' + '\t%g'*len(dat[0,:]) + '\n'
	file_max.write(fmt % ((gam,)+tuple(max_eps_cdf)))
	file_min.write(fmt % ((gam,)+tuple(min_eps_cdf)))
	file_mean.write(fmt % ((gam,)+tuple(mean_eps_cdf)))
	file_geom.write(fmt % ((gam,)+tuple(geom_eps_cdf)))
file_max.close()
file_min.close()
file_mean.close()
file_geom.close()
