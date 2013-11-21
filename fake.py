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
	return s, p1*p2*npp, p1*f2*npf, p2*f1*nfp, f1*f2*nff



## fm, pm = 0.0808 , 0.917
## fe, pe = 0.1    , 0.845
## fm, pm = 0.06 , 0.92
## fe, pe = 0.076    , 0.85
fm, pm = 0.016, 0.85
fe, pe = 0.016, 0.85

## fm, pm = 0.0409, 0.804
## fe, pe = 0.07  , 0.751
## mmntt, mmntl       , mmnll = 65 ,    696 ,    310
## emntt, emntl, emnlt, emnll = 101,    546 ,    639 ,    643
## eentt, eentl       , eenll = 75 ,    441 ,    482

## preselection
mmntt, mmntl       , mmnll = 0.81,   157.52+39.20    ,   40.40
emntt, emntl, emnlt, emnll = 6.38,   266.54,   211.77,  101.22
eentt, eentl       , eenll = 4.82,   171.60+78.56    ,   68.55 

## ## selection
## mmntt, mmntl       , mmnll = 0.000, 66.836         , 12.657
## emntt, emntl, emnlt, emnll = 1.998, 107.915, 71.055, 31.975
## eentt, eentl       , eenll = 1.110, 104.362, 27.756

mmpred = pred  (pm, fm        , mmntt, mmntl       , mmnll)
empred = predOF(pm, fm, pe, fe, emntt, emntl, emnlt, emnll)
eepred = pred  (pe, fe        , eentt, eentl       , eenll)

mmfakes = mmpred[0]
emfakes = empred[0]
eefakes = eepred[0]

mmnpp = mmpred[1]
mmnpf = mmpred[2]
mmnff = mmpred[3]


emnpp = empred[1]
emnpf = empred[2]
emnfp = empred[3]
emnff = empred[4]


eenpp = eepred[1]
eenpf = eepred[2]
eenff = eepred[3]

tot = mmfakes + emfakes + eefakes

## print 'Predicted npp, npf and nff'
## print '%4s | %4s | %4s   ||  %4s | %4s | %4s  ||  %4s | %4s | %4s' %('npp' , 'npf', 'nff', 'npp' , 'npf', 'nff', 'npp' , 'npf', 'nff')
## print '%.2f | %.2f | %.2f   ||  %.2f | %.2f | %.2f | %.2f  ||  %.2f | %.2f | %.2f' %(mmnpp, mmnpf, mmnff, emnpp, emnpf, emnfp, emnff, eenpp, eenpf, eenff)
print 
print 'mmntt: %6.2f mmntl: %6.2f mmnll: %6.2f  |  emntt: %6.2f emntl: %6.2f emnlt: %6.2f emnll: %6.2f | eentt: %6.2f eentl: %6.2f eenll: %6.2f |'            %(
mmntt, mmntl, mmnll, emntt, emntl, emnlt, emnll,  eentt, eentl, eenll)
## print 'mmnpp: %6.2f mmnpf: %6.2f mmnff: %6.2f  |  emnpp: %6.2f emnpf: %6.2f emnfp: %6.2f emnff: %6.2f | eenpp: %6.2f eenpf: %6.2f eenff: %6.2f | unweighted' %(
## mmnpp/(pm*pm), mmnpf/(fm*pm), mmnff/(fm*fm), emnpp/(pe*pm), emnpf/(pm*fe), emnfp/(pe*fm), emnff/(fm*fe),  eenpp/(pe*pe), eenpf/(pe*fe), eenff/(fe*fe))
print 'mmnpp: %6.2f mmnpf: %6.2f mmnff: %6.2f  |  emnpp: %6.2f emnpf: %6.2f emnfp: %6.2f emnff: %6.2f | eenpp: %6.2f eenpf: %6.2f eenff: %6.2f | weighted'   %(
       mmnpp,       mmnpf,       mmnff,           emnpp,       emnpf,       emnfp,       emnff,         eenpp,       eenpf,       eenff)

print ''
print ''
print 'Predicted fakes:'
print '----------------'
print '%4s | %4s | %4s' %('mumu' , 'emu'  , 'ee')
print '%4.2f | %4.2f | %4.2f' %(mmfakes, emfakes, eefakes)
print '----------------'
print 'Sum: %4.2f' %(tot)
print '----------------'
print 'Ntt for all channels:'
print '%4.2f | %4.2f | %4.2f' %(mmntt, emntt, eentt)
print 'Sum of Ntt: %.2f' %(mmntt+emntt+eentt)
