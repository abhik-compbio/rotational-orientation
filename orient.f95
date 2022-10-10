! Orientational correlation function for interfacial water molecules
! Water molecules near 5 angstorm from protein (CA atom) are considerd only
! Written by Abhik Ghosh Moulick, SNBNCBS

integer,parameter:: atomno = 369   ! change,check last atomno of complex from pdb file
integer,parameter:: solno =  12749 ! change, waternumber
integer,parameter:: maxframe =  1005 ! change, maximum frame to consider	
real xc(maxframe,1:atomno,3),xo(maxframe,1:solno,3)
real w1(maxframe,1:solno,3)
real x1(3),y1(3),z1(3)
integer ichain,atnum,gap,p,resnumb,r
character line*80
character atnam*4, junk*1
integer resnum,j,i,steps,nbin,resinum,resn,pp
real num(1:solno),distance(0:1005),counter,dist1
real dist, dist_x, dist_y, dist_z, dist_gap
real angle, sqrtx1, sqrty1
data delt/0.1/
open(10, file='../../../md_50_51ns_protein-nojump-orient.pdb') ! Change file name
open(20, file='orient-water-around-protein-free-block1.dat')! Change file name

!open(22, file='output-test.dat')
resnumb = 0
! Read pdb file
ichain=1
do while (ichain.ge.1)
read (10,'(a)') line
if (line(1:4).eq.'ATOM') then
read (line,1) resnum,atnam,junk,resn,x,y,z
!write (22,1) resnum,atnam,junk,resn,x,y,z
1 format (6x,i5,2x,a4,3x,a2,1x,i4,3x,3f8.3)
if(atnam.eq.'C' .or. atnam .eq. 'CA' .or. atnam.eq.'N' ) then
	xc(ichain,resnum,1) = x
	xc(ichain,resnum,2) = y
	xc(ichain,resnum,3) = z
	!write(*,*) atnam,resnum, xc(ichain,resnum,1)
end if

if(atnam.eq.'OW' .and. resnum .gt. atomno) then
	resnumb = resnumb + 1
	xo(ichain,resnumb,1) = x
	xo(ichain,resnumb,2) = y
	xo(ichain,resnumb,3) = z
	else if(atnam.eq.'HW1') then
	w1(ichain,resnumb,1) = x
	w1(ichain,resnumb,2) = y
	w1(ichain,resnumb,3) = z 
end if				
end if	
if (line(1:4).eq.'TER ') THEN
ichain = ichain + 1
resnumb = 0
maxchain = ichain
!write(*,*) maxchain
if (maxchain .eq. 1001) exit
end if
end do
 close(10)

! Orientational correlation function calculation
      	
! Initialize array to control interfacial water number
! Any water molecule which is at 5 angstorm distance from C-alpha atom
! of any residue should not consider for further calculations

do p = 1,solno! change loop indices
	num(p) = 0.0
end do
      	
!bin_width = 0.1  ! Define bin width

do gap = 0, 50
!write(*,*) gap
M = 950
!M = steps - gap  ! M signifies number of initial points
distance(gap) = 0.0
!write(*,*) gap
do ichain = 1,M
	dist1 = 0.0
	counter = 0.0
      	do i = 1,atomno
      		do j = 1,solno  ! change loop indices
      			dist_x = xc(ichain,i,1)-xo(ichain,j,1)
      			dist_y = xc(ichain,i,2)-xo(ichain,j,2)
      			dist_z = xc(ichain,i,3)-xo(ichain,j,3)
      			dist = sqrt(dist_x**2+dist_y**2+dist_z**2)
      			if (dist .le. 4.5) then
				if (num(j) .eq. 0.0) then
      					num(j) = num(j) + 1.0
      					angle = 0.0
      					do r = 1,3
      					x1(r) = xo(ichain,j,r)-w1(ichain,j,r)
      					y1(r) = xo(ichain+gap,j,r)-w1(ichain+gap,j,r)
      					angle = angle + x1(r)*y1(r)
      					!write(*,*)j,xo(ichain,j,r),w1(ichain,j,r)
      					end do
      					sqrtx1 = sqrt(x1(1)**2+x1(2)**2+x1(3)**2)
      					sqrty1 = sqrt(y1(1)**2+y1(2)**2+y1(3)**2)
      					angle = angle/(sqrtx1*sqrty1)
      					dist1 = dist1 + angle
      					counter = counter + 1.0
      				end if
      			end if
      			
      		end do
      	end do
if (dist1 .eq. 0.0 .or. counter .eq. 0.0) then
	distance(gap) = distance(gap)
else
	distance(gap) = distance(gap)+(dist1/counter)
end if
do p = 1,solno !change loop indices
	num(p) = 0.0
end do
end do
distance(gap) = distance(gap)/(float(M))

write(20,35)float(gap),distance(gap)/distance(0)
35 	format (f10.4,'	',f16.5)
!!write(*,*)'hi'
end do

end
