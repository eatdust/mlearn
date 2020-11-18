module ioml
  use mlprec
  implicit none

  private


  public find_boundaries, save_boundaries, read_boundaries
  public cubize_paramspace
  public save_inputdata, load_inputdata
  public delete_file, livewrite, allwrite


contains

  subroutine find_boundaries(xdata,xmin,xmax)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:), intent(out) :: xmin,xmax

    integer(ip) :: ndim, i
    
    ndim = size(xdata,1)
    
    if ((size(xmin,1).ne.ndim).or.(size(xmax,1).ne.ndim)) then
       stop 'posterior_boundaries: array mismatch!'
    endif

    do i=1,ndim
       xmin(i) = minval(xdata(i,:))
       xmax(i) = maxval(xdata(i,:))
    enddo

  end subroutine find_boundaries


  

  subroutine save_boundaries(filename,xmin,xmax,fmin,fmax)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:), intent(in) :: xmin,xmax
    real(fp), intent(in) :: fmin,fmax

    integer(ip) :: nunit
    integer(ip) :: j, ndim

    ndim = size(xmin,1)
    if (ndim.ne.size(xmax,1)) stop 'save_boundaries: xmin/xmax dim!'

    nunit = 410

    open(nunit,file=filename,status='unknown')
    write(nunit,*) ndim
    write(nunit,*)(xmin(j),j=1,ndim)
    write(nunit,*)(xmax(j),j=1,ndim)
    write(nunit,*)fmin
    write(nunit,*)fmax
    close(nunit)

  end subroutine save_boundaries



  subroutine read_boundaries(filename,xmin,xmax,fmin,fmax)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:), pointer, intent(inout) :: xmin,xmax
    real(fp), intent(inout) :: fmin,fmax
    integer(ip) :: nunit
    integer(ip) :: j, ndim

    if (associated(xmin).or.associated(xmax)) stop 'read_boundaries: already done!'

    nunit = 411

    open(nunit,file=filename,status='old')
    read(nunit,*) ndim

    allocate(xmin(ndim),xmax(ndim))

    read(nunit,*)(xmin(j),j=1,ndim)
    read(nunit,*)(xmax(j),j=1,ndim)

    read(nunit,*)fmin
    read(nunit,*)fmax


    close(nunit)

  end subroutine read_boundaries


  
  subroutine cubize_paramspace(xdata,xcubes)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:,:), intent(out) :: xcubes
    
    integer(ip) :: ndim, ndata, i
    real(fp), dimension(:), allocatable :: xmin,xmax

    ndim = size(xdata,1)
    ndata = size(xdata,2)

    if ((size(xcubes,1).ne.ndim).or.(size(xcubes,2).ne.ndata)) then
       stop 'cube_rbfparamspace: mismatch arrays!'
    endif

    allocate(xmin(ndim), xmax(ndim))

    call find_boundaries(xdata,xmin,xmax)

    do i=1,ndim       
       xcubes(i,:) = (xdata(i,:) - xmin(i))/(xmax(i)-xmin(i))
    enddo

    deallocate(xmin,xmax)

  end subroutine cubize_paramspace



  
  subroutine save_inputdata(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), intent(in) :: f
    real(fp), dimension(:,:), intent(in) :: x

    integer(ip), parameter :: nunit = 211

    integer(ip) :: i,j, ndim, ndata

    ndim = size(x,1)
    ndata = size(x,2)
    
    if (size(f,1).ne.ndata) stop 'save_posterior: size mismatch!'

    open(unit=nunit,file=filename,status='unknown')
    
    write(nunit,*) ndim, ndata
    do i=1,ndata
       write(nunit,*) f(i),(x(j,i),j=1,ndim)
    enddo
    
    close(nunit)
  end subroutine save_inputdata


  
  subroutine load_inputdata(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), pointer :: f
    real(fp), dimension(:,:), pointer :: x

    integer(ip), parameter :: nunit = 211

    integer(ip) :: i,j, ndim, ndata


    if (associated(f).or.associated(x)) then
       stop 'load_posterior: data already loaded!'
    endif

    open(unit=nunit,file=filename,status='old')
    
    read(nunit,*) ndim, ndata

    allocate(f(ndata),x(ndim,ndata))

    do i=1,ndata
       read(nunit,*) f(i),(x(j,i),j=1,ndim)
    enddo
    
    close(nunit)

  end subroutine load_inputdata



   subroutine delete_file(name)
    implicit none
    character(len=*) :: name
    logical :: isthere

    inquire(file=name,exist=isthere)

    if (isthere) then
       open(unit=10,file=name)
       close(unit=10,status='delete')
    endif

  end subroutine delete_file




  subroutine livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: istat,nunit
    real(fp) :: x,a
    real(fp), optional :: b,c,d,e,f,g
    
    nunit = 412

    open(nunit,file=name,position='append',status='unknown')

    if (.not.present(b)) then
       write(nunit,100) x,a         
    elseif (.not.present(c)) then           
       write(nunit,100) x,a,b         
    elseif (.not.present(d)) then             
       write(nunit,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(nunit,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(nunit,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(nunit,100) x,a,b,c,d,e,f            
    else
       write(nunit,100) x,a,b,c,d,e,f,g            
    endif

    close(nunit)

100 format(8(ES25.16E3))      

  end subroutine livewrite

  subroutine allwrite(name,x,a,b,c,d,e,f,g)
      implicit none
      character(*) :: name
      integer(ip) :: j,npts
      real(fp) :: x(:),a(:)
      real(fp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,position='append',status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES25.16))      

    end subroutine allwrite


  

  end module ioml
