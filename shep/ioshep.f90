module ioshep
  use prec
  implicit none

  private
  
  public save_shepdata, load_shepdata

contains


  subroutine save_shepdata(filename,rmax)
    use qshepmdata, only : A, IW
    use qshepmdata, only : DX, XMIN, RSQ, WS
    use qshepmdata, only : LCELL, LNEXT

    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(in) :: rmax
    integer :: nunit

    integer :: i,j
    integer :: n, ntt, m, nlcell

    if (.not.allocated(A)) stop 'save_shepdata: not found!'

    n = size(A,1)
    ntt = size(A,2)
    m = size(xmin)
    nlcell = size(LCELL,1)

    nunit = 511
    open(nunit,file=filename,status='unknown')
    write(nunit,*) n,ntt,m,nlcell
    write(nunit,*) rmax

    do i=1,m
       write(nunit,*) DX(i),XMIN(i)
       write(nunit,*) (IW(i,j),j=1,5)
    enddo

    do i=1,n
       write(nunit,*) RSQ(i),LNEXT(i)
       write(nunit,*) (A(i,j),j=1,ntt)
    enddo

    do i=1,ntt*ntt
       write(nunit,*) WS(i)
    enddo

    do i=1,nlcell
       write(nunit,*) LCELL(i)
    enddo

    close(nunit)

  end subroutine save_shepdata

 

  subroutine load_shepdata(filename,rmax)
    use qshepmdata, only : A, IW
    use qshepmdata, only : DX, XMIN, RSQ, WS
    use qshepmdata, only : LCELL, LNEXT

    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(out) :: rmax
    integer :: nunit

    integer :: i,j
    integer :: n, ntt, m, nlcell

    if (allocated(A)) stop 'load_shepdata: already loaded'

    n = size(A,1)
    ntt = size(A,2)
    m = size(xmin)
    nlcell = size(LCELL,1)

    nunit = 511
    open(nunit,file=filename,status='old')
    read(nunit,*) n,ntt,m,nlcell
    read(nunit,*) rmax

    allocate(A(n,ntt))
    allocate(DX(m),XMIN(m),IW(m,5))
    allocate(RSQ(N),LNEXT(N))
    allocate(LCELL(nlcell))
    allocate(WS(ntt*ntt))
    

    do i=1,m
       read(nunit,*) DX(i),XMIN(i)
       read(nunit,*) (IW(i,j),j=1,5)
    enddo

    do i=1,n
       read(nunit,*) RSQ(i),LNEXT(i)
       read(nunit,*) (A(i,j),j=1,ntt)
    enddo

    do i=1,ntt*ntt
       read(nunit,*) WS(i)
    enddo

    do i=1,nlcell
       read(nunit,*) LCELL(i)
    enddo

    close(nunit)

  end subroutine load_shepdata


end module ioshep
