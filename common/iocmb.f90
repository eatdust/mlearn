module iocmb
  use mlprec
  implicit none

  private

  integer(ip), parameter :: nrecmax = 1000000
  integer, parameter :: ndimmax = 5

  public read_binned_posterior, posterior_boundaries
  public save_posterior, load_posterior


  interface read_binned_posterior
     module procedure read_binned_posterior_vector, read_binned_posterior_matrix
  end interface read_binned_posterior


  interface load_posterior
     module procedure load_posterior_vector, load_posterior_matrix
  end interface load_posterior


  interface save_posterior
     module procedure save_posterior_vector, save_posterior_matrix
  end interface save_posterior


contains



  subroutine read_binned_posterior_matrix(filename,posterior,params,logize)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:,:),  pointer :: posterior
    real(fp), dimension(:,:), pointer :: params
    logical, intent(in), optional :: logize

    real(fp), dimension(:), pointer :: postvector

    call read_binned_posterior_vector(filename,postvector,params,logize)

    posterior(1:1,1:size(postvector,1)) => postvector(1:size(postvector,1))
    
  end subroutine read_binned_posterior_matrix



  subroutine read_binned_posterior_vector(filename,posterior,params,logize)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:),  pointer :: posterior
    real(fp), dimension(:,:), pointer :: params
    logical, intent(in), optional :: logize

    logical, save :: logpost = .true.
    integer, parameter :: nzeroskip = 0

    real(fp), save :: eps = exp(-10._fp)
    real(fp), save :: minNonZero = huge(1._fp)


    integer, parameter :: nunit = 210


    real(fp), save :: fpbuffer
    character(LEN=200) :: cbuffer
    
    real(fp), dimension(ndimmax), save :: statbuffer

    integer(ip) :: i,j,ioerr,nrec,ndim,nnzero,count

    if (present(logize)) logpost = logize


    if (associated(posterior).or.associated(params)) then
       stop 'read_binned_posterior: data already loaded!'
    endif

    write(*,*)'read_binned_posterior:'         
    write(*,*)'opening ',filename
   
!counting columns
    ndim = count_column(filename,'E') - 1
    
    
!counting number records and number of non-zero posterior values
    open(unit=nunit,file=filename,status='old')
    nnzero=0
    nrec=0
    count = 0
    do i=1,nrecmax
       read(nunit,*,iostat=ioerr) statbuffer(1:ndim+1)
       if (ioerr.ne.0) exit
       nrec = nrec + 1

       if (statbuffer(1).eq.0._fp) then
          count = count + 1
          if (count.le.nzeroskip) cycle
          count = 0
       else
          minNonZero = min(minNonZero,statbuffer(1))
       endif

       nnzero = nnzero + 1
    enddo
    close(nunit)

    write(*,*)'number of params:',ndim
    write(*,*)'number or records:',nrec
    write(*,*)'number of bins kept:',nnzero
    write(*,*)'ln(minNonZero)=    ',log(minNonZero)

!rbf
!    eps = exp(1._fp+real(int(log(minNonZero)),fp))

!shep
    eps = exp(0._fp+real(nint(log(minNonZero)),fp))


!reading non-zero records
    allocate(posterior(nnzero))
    allocate(params(ndim,nnzero))
         
    open(unit=nunit,file=filename,status='old')
    j=0
    count = 0
    do i=1,nrec
       read(nunit,*,iostat=ioerr) statbuffer(1:ndim+1)

       if (statbuffer(1).eq.0._fp) then
          count = count + 1
          if (count.le.nzeroskip) cycle
          count = 0
       end if

       j=j+1
       if (logpost) then
          posterior(j) = log(statbuffer(1)+eps)
       else
          posterior(j) = statbuffer(1)
       endif
       params(1:ndim,j) = statbuffer(2:ndim+1)
!       print *,'test',posterior(j),params(1:ndim,j)
       if (ioerr.ne.0) stop 'counting screwed!'
    enddo
    close(nunit)

  end subroutine read_binned_posterior_vector

  


  subroutine save_posterior_matrix(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:,:), intent(in) :: f
    real(fp), dimension(:,:), intent(in) :: x

    call save_posterior_vector(filename,f(1,:),x)

  end subroutine save_posterior_matrix





  subroutine save_posterior_vector(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), intent(in) :: f
    real(fp), dimension(:,:), intent(in) :: x

    integer, parameter :: nunit = 211

    integer :: i,j, ndim, ndata

    ndim = size(x,1)
    ndata = size(x,2)
    
    if (size(f,1).ne.ndata) stop 'save_posterior: size mismatch!'

    open(unit=nunit,file=filename,status='unknown')
    
    write(nunit,*) ndim, ndata
    do i=1,ndata
       write(nunit,*) f(i),(x(j,i),j=1,ndim)
    enddo
    
    close(nunit)
  end subroutine save_posterior_vector



  subroutine load_posterior_matrix(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:,:), pointer :: f
    real(fp), dimension(:,:), pointer :: x

    real(fp), dimension(:), pointer :: fvec

    call load_posterior_vector(filename,fvec,x)

    f(1:1,1:size(fvec,1)) => fvec(1:size(fvec,1))

  end subroutine load_posterior_matrix


  

  subroutine load_posterior_vector(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), pointer :: f
    real(fp), dimension(:,:), pointer :: x

    integer, parameter :: nunit = 211

    integer :: i,j, ndim, ndata


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

  end subroutine load_posterior_vector





  subroutine posterior_boundaries(xdata,xmin,xmax)
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

  end subroutine posterior_boundaries


  function count_column(filename,delimiter)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: delimiter

    logical, parameter :: display = .false.

    integer(ip) :: count_column

    integer(ip) :: stat, j, num_delim
    character :: single_byte, CR, LF, column_delimiter

    integer(ip), parameter :: nunit = 200

    LF = char(10) ! Line Feed
    CR = char(13) ! Carriage Return

    column_delimiter = delimiter    

    open (nunit, file=filename, form='unformatted', &
         access='direct', status='old', recl = 1, iostat=stat)
    if (stat .ne. 0) stop 'read_header: missing file!'

    ! process header line of the file
    j = 0
    num_delim = 0
    single_byte='o'

    do while ((single_byte .ne. CR) .and. (single_byte .ne. LF))
       j = j + 1
       read(nunit, rec = j) single_byte
       if (single_byte .eq. column_delimiter) then
          num_delim = num_delim + 1
       end if
       !write (*,'(I3,5x,I3,5x,A)') j, ichar(single_byte), single_byte
    end do
    close(nunit)

    if (display) then
       write (*,*)'delimiter ',delimiter,' found ',num_delim,'times'
    endif

    count_column = num_delim    
  end function count_column





end module iocmb
