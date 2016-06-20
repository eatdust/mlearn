program fnnmain
  use fnnprec
  use fnn
  use ioml

  implicit none

  real(fp), dimension(:,:), pointer :: xdata => null()
  real(fp), dimension(:,:), pointer :: fdata => null()
  real(fp), dimension(:,:), pointer:: xcubes => null()


  real(fp), dimension(:), pointer :: xnmin, xnmax
  real(fp), dimension(:), allocatable :: x   

  integer, parameter :: lenmax = 32
  character(len=lenmax), dimension(:), allocatable :: names

  real(fp) :: f, fmin, fmax


  integer :: i,j
  integer, parameter :: ndump = 100
  real(fp) :: lnA, sr1, sr2, sr3
  real(fp) :: lnAmin, lnAmax, sr1min, sr1max, sr2min, sr2max
  real(fp) :: sr3min,sr3max


  integer(ip), parameter :: nout = 1
  real(fp), dimension(nout) :: fnorm
  integer(ip) :: ndata, ndim


  integer(ip), parameter :: numLayers = 3
  integer(ip), dimension(numLayers) :: nNeurons
  integer(ip) :: maxEpochs, epochBetweenReports
  real(fp) :: maxError

  integer(ip) :: maxNeurons, neuronsBetweenReports

  character(len=*), parameter :: training = 'normal'


  call read_binned_posterior('test_posterior.dat',fdata,xdata)

  ndata = size(fdata,2)
  ndim = size(xdata,1)
  print *,'ndata= ',ndata
  print *,'ndim= ',ndim

  if (size(xdata,2).ne.ndata) stop 'internal error'
 
  fmax = maxval(fdata)
  fmin = minval(fdata)
  print *,'fdata max',fmax
  print *,'fdata min',fmin

  allocate(xnmin(ndim), xnmax(ndim))
  allocate(xcubes(ndim,ndata))


  call posterior_boundaries(xdata,xnmin,xnmax)

  print *,'xnmin',xnmin
  print *,'xnmax',xnmax
  print *

  call cubize_paramspace(xdata,xcubes)




  nNeurons(1) = ndim
  nNeurons(2) = 20
  nNeurons(3) = nout

  
  maxEpochs = 10000
  epochBetweenReports = 100
  maxError = 1e-3

  select case (training)

  case ('normal')
   
     call fnn_create_ann(numLayers,nNeurons)
     call fnn_set_activation_function_hidden('FANN_SIGMOID')
     call fnn_set_activation_function_output('FANN_SIGMOID')
     call fnn_print_connections()

     call fnn_create_train(ndata,ndim,nout,xcubes,fdata)
     call fnn_set_training_algorithm('FANN_TRAIN_RPROP')

     call fnn_set_train_error_function('FANN_ERRORFUNC_LINEAR')
     print *,'error function: ',fnn_get_train_error_function()

!     call fnn_set_bit_fail_limit(0.5_fp)
!     print *,'get bit', fnn_get_bit_fail_limit()
     call fnn_set_train_stop_function('FANN_STOPFUNC_BIT')


!     call fnn_set_learning_rate(0.9_fp)
     print *,'learning rate= ',fnn_get_learning_rate()

     call fnn_scale_output_train_data(0._fp,1._fp)
     
     call fnn_train_on_data(maxEpochs,epochBetweenReports,maxError)
!     call fnn_save_train('fnntrain.dat')

     call fnn_print_connections()

     print *,'mse= ', fnn_test_data()

     call fnn_save_ann('fnndata.dat')

     call save_boundaries('bounds.dat',xnmin,xnmax,fmin,fmax)


  case ('cascade')

     call fnn_create_shortcut_array(2_ip,(/nNeurons(1),nNeurons(3)/))


     call fnn_print_connections()


     call fnn_create_train(ndata,ndim,nout,xcubes,fdata)
     call fnn_scale_output_train_data(0._fp,1._fp)



     call fnn_set_train_error_function('FANN_ERRORFUNC_LINEAR')
     print *,'error function: ',fnn_get_train_error_function()

     maxNeurons = 30
     neuronsBetweenReports = 1
     maxError = 1e-3
     allocate(names(2))
     names(1) = 'FANN_SIGMOID'
     names(2) = 'FANN_GAUSSIAN'
!     names(3) = 'FANN_SIGMOID_SYMMETRIC'
!     names(4) = 'FANN_GAUSSIAN_SYMMETRIC'
     call fnn_set_cascade_activation_functions(names)
     deallocate(names)
     print *,'activation functions= ',fnn_get_cascade_activation_functions()

!     call fnn_set_cascade_max_cand_epochs(1000_ip)
!     print *,'max candidate epochs= ',fnn_get_cascade_max_cand_epochs()

     call fnn_set_cascade_activation_steepnesses((/0.1_fp,0.2_fp,0.3_fp,0.4_fp,0.5_fp &
          ,0.6_fp,0.7_fp,0.8_fp,0.9_fp/))

     call fnn_set_train_stop_function('FANN_STOPFUNC_BIT')

!     call fnn_set_bit_fail_limit(0.3_fp)
!     print *,'get bit', fnn_get_bit_fail_limit()

    call fnn_set_training_algorithm('FANN_TRAIN_RPROP')


     call fnn_cascadetrain_on_data(maxNeurons,neuronsBetweenReports,maxError)
     call fnn_print_connections()
     
     print *,'mse= ', fnn_test_data()
     call fnn_save_ann('fnndata.dat')
     call save_boundaries('bounds.dat',xnmin,xnmax,fmin,fmax)

  case default

     deallocate(fdata,xcubes)

     call fnn_create_ann('fnndata.dat')

     call fnn_print_connections()

     
  end select


!test the running  

  
  allocate(x(ndim))
  
  lnAmin = xnmin(1)
  lnAmax = xnmax(1)
  sr1min = xnmin(2)
  sr1max = xnmax(2)
  sr2min = xnmin(3)
  sr2max = xnmax(3)
  sr3min = xnmin(4)
  sr3max = xnmax(4)

  call delete_file('output.dat')
  call delete_file('output2.dat')
  call delete_file('output3.dat')
  lnA = 3.09
  sr3 = 0.0


  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)

     do j=1,ndump
        sr2 = sr2min + (sr2max-sr2min)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        fnorm = fnn_run(x)
        f = fmin + fnorm(1)*(fmax-fmin)

        call livewrite('output.dat',sr1,sr2,f)

     enddo
  enddo

  sr3=0
  sr2 = 0.04
  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)
     
     do j=1,ndump
        lnA = lnAmin + (lnAmax-lnAmin)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        fnorm = fnn_run(x)
        f = fmin + fnorm(1)*(fmax-fmin)

        call livewrite('output2.dat',sr1,lnA,f)

     enddo

  enddo

  if (ndim.lt.4) stop

  lnA = 3.09
  sr2 = 0.04
  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)
     
     do j=1,ndump
        sr3 = sr3min + (sr3max-sr3min)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        fnorm = fnn_run(x)
        f = fmin + fnorm(1)*(fmax-fmin)

        call livewrite('output3.dat',sr1,sr3,f)

     enddo

  enddo

  deallocate(xnmin,xnmax)
  deallocate(x)



  call fnn_free_ann()
  

end program fnnmain
