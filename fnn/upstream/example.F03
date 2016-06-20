program example
!also import iso_c_bindings
  use fann
  implicit none

  type(C_PTR) :: ann
  type(C_PTR) :: train

  integer, parameter :: sp = C_FLOAT
  integer, parameter :: ft = FANN_TYPE

  integer, parameter :: num_layer = 3
  integer, parameter :: nin = 3
  integer, parameter :: nout= 1
  integer, dimension(num_layer) :: layers

  integer, parameter :: ndata = 100
  real, dimension(nin,ndata) :: inTrainData
  real, dimension(nout,ndata) :: outTrainData

  integer :: max_epochs, epochs_between_reports
  real(sp) ::  desired_error

  real(ft), dimension(nin) :: x

!input
  layers(1) = nin
!hidden
  layers(2) = 20
!outout
  layers(3) = nout

!the net, with SIGMOID
  ann = fann_create_standard_array(num_layer,layers)
  call fann_set_activation_function_hidden(ann,enum_activation_function('FANN_SIGMOID'))
  call fann_set_activation_function_output(ann,enum_activation_function('FANN_SIGMOID'))

  call fann_print_connections(ann)

!my training data. Let's learn make a neural net which is a random generator :)
  call random_number(inTrainData)
  call random_number(outTrainData)
  
  train = fann_create_train_from_callback(ndata,nin,nout,C_FUNLOC(mytrain_callback))

!training
  call fann_set_training_algorithm(ann,enum_training_algorithm('FANN_TRAIN_RPROP'))

  max_epochs = 10000
  epochs_between_reports = 1000
  desired_error = 0.001
  call fann_train_on_data(ann,train,max_epochs,epochs_between_reports,desired_error)
  call fann_print_connections(ann)

!testing
  
  x = (/0.1_ft,0.5_ft,1._ft/)

!running
  print *, 'ann(x)= ',f_fann_run(ann,x)
  
!saving
  print *,'saving...', fann_save(ann,f_c_string('arg.dat'))
  call fann_destroy(ann)

!loading
  print *,'loading...'
  ann = fann_create_from_file(f_c_string('arg.dat'))
  print *, 'loaded ann(x)= ',f_fann_run(ann,x)

contains


 subroutine mytrain_callback(num, num_input, num_output, input, output) bind(C)
    implicit none

    integer(C_INT), value :: num, num_input, num_output
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(0:num_input-1) :: input
    integer(FANN_TYPE), dimension(0:num_output-1) :: output

    input(0:num_input-1) = int(inTrainData(1:num_input,num+1),FANN_TYPE)
    output(0:num_output-1) = int(outTrainData(1:num_output,num+1),FANN_TYPE)

#else
    real(FANN_TYPE), dimension(0:num_input-1) :: input
    real(FANN_TYPE), dimension(0:num_output-1) :: output

    input(0:num_input-1) = real(inTrainData(1:num_input,num+1),FANN_TYPE)
    output(0:num_output-1) = real(outTrainData(1:num_output,num+1),FANN_TYPE)

#endif
    
  end subroutine  mytrain_callback



end program example
