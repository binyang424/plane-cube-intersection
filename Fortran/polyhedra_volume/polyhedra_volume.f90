      program polyhedra_volume
      implicit none

! computes the volume of a polyhedron

      character*80       :: string,infile
      logical            :: lexist
      logical, parameter :: idebug = .false.
      integer            :: i,j
      integer, parameter :: nmax = 200
      integer            :: nvert,nfaces,face(nmax),nface_vertices(nmax),ivertices(nmax,nmax)
      real*8             :: x(nmax),y(nmax),z(nmax), &
                            xp(nmax),yp(nmax),zp(nmax), &
                            px, py, pz, xn, &
                            volume,dv, &
                            manswer,diff

! popular formats
09    format(a)
10    format(a,1p3e14.6)
11    format(a,i3,a,1pe9.2,a,1pe9.2,a,1pe9.2)


! get the file name from the command line if it was passed
      infile = ' '
      call getarg(1, infile)

      if (infile .eq. ' ') then
       write(6,*) 'give file ndate with vertex and face info =>'
       read(5,09) infile
      end if

      inquire (file=infile,exist=lexist)
      if (.not.lexist) stop 'file does not exist'



! open and read header
      open(unit=2,file=trim(infile),status='old')
      read(2,09) string          ! header info
      read(2,*) manswer


! number of vertices and their 3d coordinates
      read(2,*) nvert
      if (nvert .gt. nmax) stop 'nvert > nmax'
      read(2,*) (x(i),y(i),z(i), i=1,nvert)

      if (idebug) write(6,'(i4,3f9.2)') (i,x(i),y(i),z(i), i=1,nvert)

! number of faces, face number and list of vertices associated with that face number
      read(2,*) nfaces
      if (nfaces .gt. nmax) stop 'nfaces > nmax'

      if (idebug) write(6,*) 'nfaces ', nfaces

      do i=1,nfaces
       read(2,*) face(i),nface_vertices(i)
       read(2,*) (ivertices(j,i), j=1,nface_vertices(i))
       if (idebug) write(6,'(99i3)') face(i),nface_vertices(i), (i,ivertices(j,i), j=1,nface_vertices(i))
      enddo
      close(unit=2)


! a few sanity checks on the face data
      do i=1,nfaces
       do j=1,nface_vertices(i)
        if (ivertices(j,i) .gt. nvert) then
         write(6,*) 
         write(6,*) 'nfaces ', i, ' vertex ',j,' has a value', ivertices(j,i), ' larger than the number of vertices ', nvert
         stop 'bad face listing'
        end if
        if (ivertices(j,i) .lt. 1) then
         write(6,*) 
         write(6,*) 'nfaces ', i, ' vertex ',j,' has a value', ivertices(j,i), ' is less than 1'
         stop 'bad face listing'
        end if
       end do
      enddo



! the volume should be independent of the common "height" used for each polygon
!
! choosing any point on the surface as the height, say the first point, 
!            px = x(1)  ; py = y(1) ; pz = z(1)
! works fine for convex shapes, but fails for non-convex shapes such as stellated polyhedra.
!
! choosing an interior point as the height, say the "centroid", works fine for convex and non-convex polyhedra. 
! if only vertices have weight, then the mean value of the weighted sum of all the points is the "centroid: 
! r_c = sum (r_i * m_i) / nvert. for unity mass weightings m, the even simpler r_c = sum (r_i) / nvert. 
! we'll use this interior point as the height. this effectly selects the chosen origin point, an interior point, for symmetric polyhedra.

      xn = 1.0d0/float(nvert)
      px = sum(x(1:nvert)) * xn
      py = sum(y(1:nvert)) * xn
      pz = sum(z(1:nvert)) * xn

      write(6,*) 
      write(6,09) trim(infile)

      volume = 0.0d0
      do i = 1,nfaces

! create a polygon for each face
       do j=1,nface_vertices(i)
        xp(j) = x(ivertices(j,i))
        yp(j) = y(ivertices(j,i))
        zp(j) = z(ivertices(j,i))
       enddo

! get the volume associated with this polygon

       call pyramid_volume(xp,yp,zp,nface_vertices(i),px,py,pz,dv)

       volume = volume + dv
      enddo

!      write(6,*)
      write(6,12) 'volume = ',volume,'   mathematica ',manswer,'   difference ',volume-manswer
12    format(a,f20.16,a,f20.16,a,1pe24.16)


      end program polyhedra_volume



      include 'polylib.f90'
