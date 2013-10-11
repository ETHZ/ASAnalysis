#!/usr/bin/python

# standalone function for calculating
# the fake estimates


def pred(p, f, ntt, ntl, nll):
	c   = 1/((p-f)*(p-f))
	npp = c * ( ntt*(1-f)*(1-f)      - ntl*f*(1-f)             + nll*f*f   )
	npf = c * ( ntt*(-2)*(1-p)*(1-f) + ntl*(f*(1-p) + p*(1-f)) - nll*2*f*p )
	nff = c * ( ntt*(1-p)*(1-p)      + ntl*p*(1-p)             + nll*p*p   )
	s = p*f*npf + f*f*nff
	return s, p*p*npp, f*p*npf, f*f*nff

	
def predOF(p1, f1, p2, f2, ntt, ntl, nlt, nll):
	npp = ((f1-p1)*(f2-p2))**(-1) * ( ntt*(f1-1)*(f2-1) + ntl*(f1-1)*f2 + nlt*f1*(f2-1) + nll*f1*f2 )
	npf = ((f1-p1)*(f2-p2))**(-1) * ( ntt*(f1-1)*(1-p2) - ntl*(f1-1)*p2 + nlt*f1*(1-p2) - nll*f1*p2 )
	nfp = ((f1-p1)*(f2-p2))**(-1) * ( ntt*(1-p1)*(f2-1) + ntl*(1-p1)*f2 - nlt*p1*(f2-1) - nll*p1*f2 )
	nff = ((f1-p1)*(f2-p2))**(-1) * ( ntt*(1-p1)*(1-p2) - ntl*(1-p1)*p2 - nlt*p1*(1-p2) + nll*p1*p2 )
	s = p1*f2*npf + f1*p2*nfp + f1*f2*nff
	return s, p1*p2*npp, p1*f1*npf, p2*f1*nfp, f1*f2*nff



## fm, pm = 0.0808 , 0.917
## fe, pe = 0.1    , 0.845
## fm, pm = 0.06 , 0.92
## fe, pe = 0.076    , 0.85
fm, pm = 0.0764 , 0.916
fe, pe = 0.0885 , 0.82

mmntt, mmntl       , mmnll = 2, 14     , 2
emntt, emntl, emnlt, emnll = 5, 51,  12, 10
eentt, eentl       , eenll = 4, 32     , 6

## mmntt, mmntl       , mmnll = 111,  1005    , 226
## emntt, emntl, emnlt, emnll = 220,  807, 852, 610
## eentt, eentl       , eenll = 146,  634     , 395

mmpred = pred  (pm, fm        , mmntt, mmntl       , mmnll)
empred = predOF(pm, fm, pe, fe, emntt, emntl, emnlt, emnll)
eepred = pred  (pe, fe        , eentt, eentl       , eenll)

mmfakes = mmpred[0]
emfakes = empred[0]
eefakes = eepred[0]

mmnpp = mmpred[1]
emnpp = empred[1]
eenpp = eepred[1]

mmnpf = mmpred[2]
emnpf = empred[2]
emnfp = empred[3]
eenpf = eepred[2]

mmnff = mmpred[3]
emnff = empred[4]
eenff = eepred[3]

tot = mmfakes + emfakes + eefakes

## print 'Predicted npp, npf and nff'
## print '%4s | %4s | %4s   ||  %4s | %4s | %4s  ||  %4s | %4s | %4s' %('npp' , 'npf', 'nff', 'npp' , 'npf', 'nff', 'npp' , 'npf', 'nff')
## print '%.2f | %.2f | %.2f   ||  %.2f | %.2f | %.2f | %.2f  ||  %.2f | %.2f | %.2f' %(mmnpp, mmnpf, mmnff, emnpp, emnpf, emnfp, emnff, eenpp, eenpf, eenff)
print 'Predicted fakes:'
print '----------------'
print '%4s | %4s | %4s' %('mumu' , 'emu'  , 'ee')
print '%4.2f | %4.2f | %4.2f' %(mmfakes, emfakes, eefakes)
print '----------------'
print 'Sum: %4.2f' %(tot)
