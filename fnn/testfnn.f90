program testfnn
  use prec, only : ip, fp
  use fnn
  implicit none

  integer(ip), parameter :: ndim = 5
  integer(ip), parameter :: nout = 1
  
  integer(ip), parameter :: numLayers = 3
  integer(ip), dimension(numLayers) :: nNeurons

  integer(ip) :: ndata

  integer(ip) :: maxNeurons, maxEpoch, epochBetweenReport, neuronsBetweenReports
  real(fp) :: maxError, error, slope

  integer(ip), dimension(:), allocatable :: ibuffer
  type(fann_connection), dimension(:), allocatable :: connections

  real(fp), dimension(:,:), allocatable :: xcubes, fdata

  real(fp), dimension(:), allocatable :: invector
  real(fp), dimension(:), allocatable :: outvector

  character(len=40), dimension(:), allocatable :: names

  nNeurons(1) = ndim
  nNeurons(2) = 5
  nNeurons(3) = nout
   
  ndata = 100

  call fnn_create_ann(numLayers,nNeurons)

  call fnn_print_connections()

!  call fnn_print_parameters()

  print *,'input output neurons= ',fnn_get_num_input(), fnn_get_num_output(), fnn_get_total_neurons()
  print *,'conn type= ',fnn_get_total_connections(), fnn_get_network_type()
  print *,'rate= layers= ',fnn_get_connection_rate(), fnn_get_num_layers()
!  print *,'decimal point= ',fnn_get_decimal_point()
!  print *,'multiplier= ',fnn_get_multiplier()

  allocate(ibuffer(3))
  call fnn_get_layer_array(ibuffer)
  print *,'layers= ',ibuffer

  call fnn_get_bias_array(ibuffer)
  print *,'biases= ',ibuffer
  deallocate(ibuffer)



  allocate(connections(fnn_get_total_connections()))

  call fnn_get_connection_array(connections)
  print *,'connections(1)', connections(1)

  connections(1)%weight = 0.12345
  call fnn_set_weight_array(connections)
  call fnn_get_connection_array(connections)
  print *,'connections(1)', connections(1)

  call fnn_set_weight(0_ip,6_ip,0.1_fp)
  call fnn_get_connection_array(connections)
  print *,'connections(1)', connections(1)

  deallocate(connections)

  allocate(xcubes(ndim,ndata),fdata(nout,ndata))
  call random_number(xcubes)
  call random_number(fdata)

  allocate(invector(ndim), outvector(nout))

  print *,'mse= bit=',fnn_get_mse(), fnn_get_bit_fail()
  call fnn_train(xcubes(:,1),fdata(:,1))
  print *,'mse= bit=',fnn_get_mse(), fnn_get_bit_fail()
  call fnn_train(xcubes(:,2),fdata(:,2))
  print *,'mse= bit=',fnn_get_mse(), fnn_get_bit_fail()
  print *,'fnn_test= ',fnn_test(xcubes(:,2),fdata(:,2)), fdata(:,2)


  call fnn_randomize_weights(-0.1_fp,0.1_fp)

  call fnn_create_train(ndata, ndim, nout)
  call fnn_free_train()

 
  call fnn_create_train(ndata,ndim,nout,xcubes,fdata)
  call fnn_shuffle_train_data()

  call fnn_set_input_scaling_params(-1._fp,1._fp)
  invector = 1._fp
  call fnn_scale_input(invector)
  print *,'invector rescale',invector
  call fnn_descale_input(invector)
  print *,'invector descale',invector

  call fnn_set_output_scaling_params(-2._fp,2._fp)
  outvector = 1._fp
  call fnn_scale_output(outvector)
  print *,'outvector rescale',outvector

  call fnn_descale_output(outvector)
  print *,'outvector descale',outvector
  

  call fnn_set_scaling_params(0._fp,1._fp,0._fp,1._fp)
  call fnn_clear_scaling_params()


  error = fnn_test_data()

  print *,'test error= ', error

  call fnn_save_train('train.dat')

  call fnn_save_train_to_fixed('train_fixed.dat',16_ip)

  call fnn_init_weights()

  call fnn_read_train('train.dat')

  call fnn_scale_input_train_data(0._fp,2._fp)
  call fnn_scale_output_train_data(0._fp,3._fp)
  call fnn_scale_train_data(0._fp,1._fp)

  print *,'length train data ',fnn_length_train_data()

  print *,'train nin= nout= ',fnn_num_input_train_data(), fnn_num_output_train_data()

  print *,'fdata= fnn_run= ',fdata(1,2), fnn_run(xcubes(:,2))

  print *,'ALGO',fnn_get_training_algorithm()

  call fnn_set_training_algorithm('FANN_TRAIN_RPROP')

  print *,'ALGO',fnn_get_training_algorithm()

  print *,'RATE',fnn_get_learning_rate()

  call fnn_set_learning_rate(0.8_fp)

  print *,'RATE',fnn_get_learning_rate()

  print *,'ACT',fnn_get_activation_function(1_ip,2_ip)

  call fnn_set_activation_function('FANN_SIGMOID',1_ip,2_ip)
  
  print *,'ACT',fnn_get_activation_function(1_ip,2_ip)

  call fnn_set_activation_function_layer('FANN_COS',1_ip)

  print *,'ACT',fnn_get_activation_function(1_ip,2_ip)

  call fnn_set_activation_function_hidden('FANN_SIN')

  print *,'ACT',fnn_get_activation_function(1_ip,2_ip)

  call fnn_set_activation_function_output('FANN_GAUSSIAN')

  print *,'ACT',fnn_get_activation_function(2_ip,1_ip)

  call fnn_set_activation_steepness(0.3_fp,1_ip,2_ip)
  slope = fnn_get_activation_steepness(1_ip,2_ip)
  print *, 'STEEP', slope 


  call fnn_set_activation_steepness_layer(0.4_fp,1_ip)
  slope = fnn_get_activation_steepness(1_ip,2_ip)
  print *, 'STEEP', slope

  call fnn_set_activation_steepness_hidden(0.5_fp)
  slope = fnn_get_activation_steepness(1_ip,2_ip)
  print *, 'STEEP', slope

  call fnn_set_activation_steepness_output(0.45_fp)
  slope = fnn_get_activation_steepness(2_ip,1_ip)
  print *, 'STEEP', slope

  print *,'ERRFUNC',fnn_get_train_error_function()
  call fnn_set_train_error_function('FANN_ERRORFUNC_LINEAR')
  print *,'ERRFUNC',fnn_get_train_error_function()

  print *,'STOPFUNC',fnn_get_train_stop_function()
  call fnn_set_train_stop_function('FANN_STOPFUNC_BIT')
  print *,'STOPFUNC',fnn_get_train_stop_function()

  print *, 'BITLIMIT',fnn_get_bit_fail_limit()
  call fnn_set_bit_fail_limit(0.2_fp)
  print *, 'BITLIMIT',fnn_get_bit_fail_limit()

  print *, 'QPDECAY', fnn_get_quickprop_decay()
  call fnn_set_quickprop_decay(-0.0002_fp)
  print *, 'QPDECAY', fnn_get_quickprop_decay()

  print *, 'QPMU', fnn_get_quickprop_mu()
  call fnn_set_quickprop_mu(1.8_fp)
  print *, 'QPMU', fnn_get_quickprop_mu()

  print *, 'RPIF',fnn_get_rprop_increase_factor()
  call fnn_set_rprop_increase_factor(1.3_fp)
  print *, 'RPIF',fnn_get_rprop_increase_factor()

  print *, 'RPDF',fnn_get_rprop_decrease_factor()
  call fnn_set_rprop_decrease_factor(0.4_fp)
  print *, 'RPDF',fnn_get_rprop_decrease_factor()

  print *, 'RPDMin',fnn_get_rprop_delta_min()
  call fnn_set_rprop_delta_min(0.01_fp)
  print *, 'RPDMin',fnn_get_rprop_delta_min()

  print *, 'RPDMax',fnn_get_rprop_delta_max()
  call fnn_set_rprop_delta_max(49._fp)
  print *, 'RPDMax',fnn_get_rprop_delta_max()

  print *, 'RPD0',fnn_get_rprop_delta_zero()
  call fnn_set_rprop_delta_zero(0.2_fp)
  print *, 'RPD0',fnn_get_rprop_delta_zero()


  print *,'SRWDS',fnn_get_sarprop_weight_decay_shift()
  call fnn_set_sarprop_weight_decay_shift(-7.644_fp)
  print *,'SRWDS',fnn_get_sarprop_weight_decay_shift()

  print *,'SRSETF', fnn_get_sarprop_step_error_threshold_factor()
  call fnn_set_sarprop_step_error_threshold_factor(0.2_fp)
  print *,'SRSETF', fnn_get_sarprop_step_error_threshold_factor()
  
  print *,'SRSES', fnn_get_sarprop_step_error_shift()
  call fnn_set_sarprop_step_error_shift(1.485_fp)
  print *,'SRSES', fnn_get_sarprop_step_error_shift()


  print *,'SRTemp', fnn_get_sarprop_temperature()
  call fnn_set_sarprop_temperature(0.025_fp)
  print *,'SRTemp', fnn_get_sarprop_temperature()

  maxEpoch = 10000
  epochBetweenReport = 1000
  maxError = 0.02

  call fnn_train_on_data(maxEpoch,epochBetweenReport,maxError)

  call fnn_free_ann()

  print *
  print *, 'testing cascade...(press enter)'
  print *

  call fnn_create_shortcut_array(2_ip,(/nNeurons(1),nNeurons(3)/))

  maxNeurons = 100
  neuronsBetweenReports = 1
  maxError = 0.05

  call fnn_cascadetrain_on_data(maxNeurons,neuronsBetweenReports,maxError)

  print *,'output change fraction= ', fnn_get_cascade_output_change_fraction()
  call fnn_set_cascade_output_change_fraction(0.1_fp)
  print *,'output change fraction= ', fnn_get_cascade_output_change_fraction()

  print *,'output stag epoch= ', fnn_get_cascade_output_stagnation_epochs()
  call fnn_set_cascade_output_stagnation_epochs(10_ip)
  print *,'output stag epoch= ', fnn_get_cascade_output_stagnation_epochs()

  print *,'cand change fraction= ',fnn_get_cascade_candidate_change_fraction()
  call fnn_set_cascade_candidate_change_fraction(0.1_fp)
  print *,'cand change fraction= ',fnn_get_cascade_candidate_change_fraction()

  print *,'cad stag epoch= ',fnn_get_cascade_candidate_stagnation_epochs()
  call fnn_set_cascade_candidate_stagnation_epochs(14_ip)
  print *,'cad stag epoch= ',fnn_get_cascade_candidate_stagnation_epochs()

  print *,'weight multiplier= ',fnn_get_cascade_weight_multiplier()
  call fnn_set_cascade_weight_multiplier(0.35_fp)
  print *,'weight multiplier= ',fnn_get_cascade_weight_multiplier()

  print *,'cand limit= ',fnn_get_cascade_candidate_limit()
  call fnn_set_cascade_candidate_limit(990._fp)
  print *,'cand limit= ',fnn_get_cascade_candidate_limit()

  print *,'max out epoch= ',fnn_get_cascade_max_out_epochs()
  call fnn_set_cascade_max_out_epochs(100_ip)
  print *,'max out epoch= ',fnn_get_cascade_max_out_epochs()

  print *,'min out epoch= ',fnn_get_cascade_min_out_epochs()
  call fnn_set_cascade_min_out_epochs(40_ip)
  print *,'min out epoch= ',fnn_get_cascade_min_out_epochs()

  print *,'max cand epoch= ',fnn_get_cascade_max_cand_epochs()
  call fnn_set_cascade_max_cand_epochs(140_ip)
  print *,'max cand epoch= ',fnn_get_cascade_max_cand_epochs()

  print *,'min cand epoch= ',fnn_get_cascade_min_cand_epochs()
  call fnn_set_cascade_min_cand_epochs(30_ip)
  print *,'min cand epoch= ',fnn_get_cascade_min_cand_epochs()

  print *,'num cand= ',fnn_get_cascade_num_candidates()
  print *,'acti func count= ',fnn_get_cascade_activation_functions_count()
  print *,'acti func are= ',fnn_get_cascade_activation_functions()

  allocate(names(6))
  names(6)='FANN_SIGMOID'
  names(5)='FANN_SIGMOID_SYMMETRIC'
  names(4)='FANN_GAUSSIAN'
  names(3)='FANN_GAUSSIAN_SYMMETRIC'
  names(2)='FANN_ELLIOT'
  names(1)='FANN_ELLIOT_SYMMETRIC'

  call fnn_set_cascade_activation_functions(names)
  print *,'acti func are= ',fnn_get_cascade_activation_functions()
  
  print *,'acti steep count= ',fnn_get_cascade_activation_steepnesses_count()
  print *,'acti steeps = ',fnn_get_cascade_activation_steepnesses()
  call fnn_set_cascade_activation_steepnesses((/0.2_fp,0.4_fp,0.8_fp/))
  print *,'acti steep count= ',fnn_get_cascade_activation_steepnesses_count()
  print *,'acti steeps = ',fnn_get_cascade_activation_steepnesses()

  print *,'acti groups= ',fnn_get_cascade_num_candidate_groups()
  call fnn_set_cascade_num_candidate_groups(4)
  print *,'acti groups= ',fnn_get_cascade_num_candidate_groups()



  deallocate(xcubes,fdata,names)

  

end program testfnn
