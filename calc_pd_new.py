from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19
gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


# cross sections for A >= 12 (TALYS)
ddir2 = 'tables/PD_Elastic/new/'
isotopes2 = genfromtxt(ddir2 + 'isotopes.txt')
x = genfromtxt(ddir2+'eps_elastic.txt') * eV * 1e6  # [J]
n = len(x)
d2elastic = genfromtxt(ddir2+'xs_elastic.txt',
    dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
eps2 = iR.romb_pad_logspaced(x, 513)  # padding
xs2elastic = array([iR.romb_pad_zero(x, 513) for x in d2elastic['xs']])*1e-31


fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()]


for field in fields:
    print field.name

    # Calculate total interaction rate
    R4 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2elastic])

    # save elastic
    fname = ddir2+'pd_elastic_new_%s.txt' % field.name
    output = r_[
        c_[d2elastic['Z'], d2elastic['N'], R4]]
    fmt = '%i\t%i' + '\t%g'*len(gamma)
    hdr = 'Photo-disintegration with the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)

    # calc cumulative rate for elastic scattering
    fname = ddir2+'pd_elastic_new_CDF_%s.txt' % field.name
    hdr = '# Photo-disintegration with the %s\n# Z, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps\n' % field.info
    f = open(fname,'w')
    f.write(hdr)
    fmt = '%i\t%i\t%g' + '\t%g'*len(eps2) + '\n'
    for i,x in enumerate(xs2elastic):
        C1 = iR.cumulative_rate_gamma_eps(eps2,x,gamma,field)
        for j in range(0,len(gamma)):
            r = C1[j,:]/max(C1[j,:])
            r[isnan(r)] = 0
            f.write(fmt % ( (d2elastic['Z'][i],d2elastic['N'][i],log10(gamma[j])) + tuple(r)))
    f.close()
