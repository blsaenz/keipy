  module kei_biocommon

    use kei_parameters

    implicit none
        
    ! bio flx
    real, save :: &
    	bioflux(nsb,nsb,nz), &
    	biosave(nsb,nsb,nz), &
    	biomix(0:nz,nsb,2), &
    	bioadvect(0:nz,nsb)
    	
    ! bio rate
    real, save :: &
    	amort(nsb), &
    	ex(nsb), &
    	ex1(nsb), &
    	resp(nsb), &
    	sresp(nsb), &
    	sink(nsb), &
    	uptake(nsb)
    	
    ! flow dest
    real, save :: &
    	eg(nsb,nsb), &
    	peg(nsb,nsb), &
    	pex(nsb,nsb), &
    	pmort(nsb,nsb)
    	
    ! feeding
    real, save :: &
			p(nsb,nsb), &
			pi(nsb,nsb)
!			kpsw(nsb) -- apparently not used
 			
    ! bio paras
    real, save :: &
    	param(nsb,nsb,2), &
    	itype(nsb,nsb), &
    	jpin(nsb,nsb), &
    	pin(nsb,nsb,2), &
    	jpin1(nsb,nsb), &
    	flim(nsb,nsb), &
			tc(nsb)
		integer, save :: &
			biotype(nsb)

    ! light par
    real, save :: &
    akw, &
    akc, &
    jclou, &
    parc, &
    alpha(nsb), &
    Epar(0:nz), &
    kpar(nz)
        
  end module kei_biocommon


