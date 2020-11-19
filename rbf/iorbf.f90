module iorbf
  use rbfprec, only : ip, fp
  implicit none

  private

  interface load_weights
     module procedure load_weights_array, load_weights_pointer
  end interface load_weights

  interface load_centres
     module procedure load_centres_array, load_centres_pointer
  end interface load_centres
  


  
  public save_weights, load_weights
  public save_centres, load_centres


contains



  subroutine save_weights(filename,scale,weight)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(in) :: scale
    real(fp), dimension(:), intent(in) :: weight
    
    integer :: nunit
    integer(ip) :: i, nw
   
    nw = size(weight,1)

    nunit = 310

    open(nunit,file=filename,status='unknown')
    write(nunit,*) nw, scale
    do i=1,nw
       write(nunit,*) weight(i)
    enddo
    close(nunit)
  end subroutine save_weights


  

  subroutine load_weights_array(filename,scale,weight)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(out) :: scale
    real(fp), dimension(:), allocatable, intent(inout) :: weight

    integer :: nunit
    integer(ip) :: i, nw

    if (allocated(weight)) then
       stop 'load_weight_array: weight already loaded!'
    endif
    
    nunit = 311

    open(nunit,file=filename,status='old')
    read(nunit,*) nw, scale

    write(*,*)'load_weights_array: nw= scale= ',nw, scale

    allocate(weight(nw))    
    do i=1,nw
       read(nunit,*) weight(i)
    enddo
    close(nunit)

  end subroutine load_weights_array


   subroutine load_weights_pointer(filename,scale,weight)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(out) :: scale
    real(fp), dimension(:), pointer, intent(inout) :: weight

    integer :: nunit
    integer(ip) :: i, nw

    if (associated(weight)) then
       stop 'load_weight_pointer: weight already loaded!'
    endif
    
    nunit = 311

    open(nunit,file=filename,status='old')
    read(nunit,*) nw, scale

    write(*,*)'load_weights_pointer: nw= scale= ',nw, scale

    allocate(weight(nw))    
    do i=1,nw
       read(nunit,*) weight(i)
    enddo
    close(nunit)

  end subroutine load_weights_pointer



  subroutine save_centres(filename,xctrs)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:,:), intent(in) :: xctrs

    integer :: nunit
    integer(ip) :: i, j, nw, ndim

    ndim = size(xctrs,1)
    nw = size(xctrs,2)

    nunit = 320

    open(nunit,file=filename,status='unknown')
    write(nunit,*) ndim, nw
    do j=1,nw
       write(nunit,*) (xctrs(i,j),i=1,ndim)
    enddo
    close(nunit)

  end subroutine save_centres



  subroutine load_centres_pointer(filename,xctrs)
    implicit none
    character(len=*), intent(in) :: filename    
    real(fp), dimension(:,:), pointer, intent(inout) :: xctrs

    integer :: nunit
    integer(ip) :: i, j,nw ,ndim

    if (associated(xctrs)) then
       stop 'load_centres: centres already loaded!'
    endif
    
    nunit = 341

    open(nunit,file=filename,status='old')
    read(nunit,*) ndim, nw
      
    write(*,*)'load_centres_pointer: ndim= nw= ',ndim,nw

    allocate(xctrs(ndim,nw))    
    do j=1,nw
       read(nunit,*) (xctrs(i,j),i=1,ndim)
    enddo
    close(nunit)

  end subroutine load_centres_pointer


  
  subroutine load_centres_array(filename,xctrs)
    implicit none
    character(len=*), intent(in) :: filename    
    real(fp), dimension(:,:), allocatable, intent(inout) :: xctrs

    integer :: nunit
    integer(ip) :: i, j,nw ,ndim

    if (allocated(xctrs)) then
       stop 'load_centres: centres already loaded!'
    endif
    
    nunit = 341

    open(nunit,file=filename,status='old')
    read(nunit,*) ndim, nw
      
    write(*,*)'load_centres_array: ndim= nw= ',ndim,nw

    allocate(xctrs(ndim,nw))    
    do j=1,nw
       read(nunit,*) (xctrs(i,j),i=1,ndim)
    enddo
    close(nunit)

  end subroutine load_centres_array
  

end module iorbf
