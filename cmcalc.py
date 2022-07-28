def cmcalc(natm, atmtyp, xc, yc, zc, cm, rcm, mx, m, rcm_on):
  """
  integer :: i,j,m,natm,mx,c
  character(len=2) :: atmtyp(natm)
  real*8 :: xc(natm),yc(natm),zc(natm),z(natm),dist
  real*8, dimension (natm,natm) :: cm
  real*8, dimension (m) :: rcm
  logical :: rcm_on
  """

  #! convert atom type list to atomic number Z
  z = atomtype2z(atomtype)
 
  #! sort by atomic number (z) then loop over sorted indices to create CM
  #! Implement this here <----------

  #! diagonal elements
  for i in range(natm):
      cm[i, i] = 0.5 * z[i] ** 2.4
  
  #! OFF DIAGONALS
  for i in range(natm):
      for j in range(i, natm):
      #do j=i+1,natm    
      dist = ( (xc[i]-xc[j])**2 + (yc[i]-yc[j])**2 + (zc[i]-zc[j])**2 )**.5
      if (dist > 1.0e-8):
        cm[i,j] = z[i]*z[j]/dist
      else:
        cm[i,j] = 0.0
      #! other side of diagonal is equivalent
        cm(j,i) = cm(i,j) 
  """ 
  #! CALCULATE REDUCED CM OR MODIFIED FULL CM
  106 format(f12.8)
  !OPEN(7,FILE='rcm.dat',FORM="FORMATTED",STATUS="replace")
  !write(7,*) 'CM: ',idx
  if (rcm_on) then
    rcm = cm(1:m,1)
  else
    c=0
    do i=1,mx-1
      do j=i+1,mx
        c=c+1
        rcm(c) = cm(i,j)
      enddo
    enddo
  endif
  rcm(1) = 0.0d0  ! first element is absorbing atom Z^2, set to 0 (for reduced or full CM)

  !do i=1,m
    !write(7,106) rcm[i]
  !enddo
  !close(7)
  """

    return cm
