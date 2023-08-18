module utilities
    implicit none
 
    contains

    function ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      !     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      !     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
      !     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
  
      if(idum.lt.0.or.iff.eq.0)then
      iff=1
      mj=MSEED-iabs(idum)
      mj=mod(mj,MBIG)
      ma(55)=mj
      mk=1
      do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
      end do
      do  k=1,4
          do  i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
      end do
      inext=0
      inextp=31
      idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
  
      return
  
   end function ran3

    
    function dot_product(a, b) result(dp)
       real(8), intent(in) :: a(3), b(3)
       real(8) :: dp
 
       dp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
    end function dot_product
 
    function cross_product_component(a, b, index) result(cpc)
       real(8), intent(in) :: a(3), b(3)
       integer, intent(in) :: index
       real(8) :: cpc
 
       select case(index)
             case(1) ! i-component
                cpc = a(2)*b(3) - a(3)*b(2)
             case(2) ! j-component
                cpc = a(3)*b(1) - a(1)*b(3)
             case(3) ! k-component
                cpc = a(1)*b(2) - a(2)*b(1)
             case default
                cpc = 0.0_8
       end select
    end function cross_product_component

       
    subroutine CompactArrays()
      ! Compact the arrays based on active particles
      use share
      implicit none
      npartActive = 0
      do j = 1, npart
         if (isActive(j)) then
            npartActive = npartActive + 1
            x(npartActive) = x(j)
            y(npartActive) = y(j)
            z(npartActive) = z(j)
            vx(npartActive) = vx(j)
            vy(npartActive) = vy(j)
            vz(npartActive) = vz(j)
         end if
      end do
      npart = npartActive
    end subroutine CompactArrays
    
 
    function sheath_field(z, alpha, larmor_over_debye_radius) result(sheath)
       real(8), intent(in) :: z, alpha, larmor_over_debye_radius
       real(8) :: sheath, potential
       real(8) :: C1, C2
 
       C1 = -0.00281407_8 * alpha - 2.31655435_8
       C2 = 0.00640402_8 * alpha + 0.01023915_8
 
       potential = C1 * exp(-C2 * z * larmor_over_debye_radius)
       sheath = C2 * larmor_over_debye_radius * potential
 
    end function sheath_field


 
 end module utilities