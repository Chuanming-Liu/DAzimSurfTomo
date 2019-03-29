c  The following parameters should be changed when the levels in the
c  model increase.
c	level: number of levels used in model parameterization
c	ndata: level+8
c	lev1:  level*2
c	lev2:  ndata*6
c	lev3:  ndata*6+8
c	lev5:  should be larger than ndata

	parameter(level=50000,ndata=level+8,lev1=2*level,lev2=6*level,
     1	  lev3=lev2+8,lev5=55000)
