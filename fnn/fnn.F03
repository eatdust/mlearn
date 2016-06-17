module fnn
  use fnnprec
  use fann

  implicit none

  private

!the artificial neural net
  type(C_PTR) :: ann

!training data
  type(C_PTR) :: traindata


!temporary arrays for buffering data training
  real(fp), allocatable, dimension(:,:) :: inTrainData, outTrainData



  interface fnn_create_ann
     module procedure fnn_create_standard_array, fnn_create_sparse_array, fnn_create_from_file
  end interface fnn_create_ann


  interface fnn_save_ann
     module procedure fnn_save
  end interface fnn_save_ann


  interface fnn_create_train
     module procedure fnn_create_train_empty, fnn_create_train_from_callback
  end interface fnn_create_train


  interface fnn_read_train
     module procedure fnn_read_train_from_file
  end interface fnn_read_train


  interface fnn_save_train
     module procedure fnn_save_train_to_fixed, fnn_save_train_to_file
  end interface fnn_save_train


  interface fnn_cascade_train
     module procedure fnn_cascadetrain_on_data, fnn_cascadetrain_on_file
  end interface fnn_cascade_train

  

  public fnn_save_ann, fnn_save_to_fixed, fnn_create_ann, fnn_run, fann_connection
  public fnn_create_standard_array, fnn_create_sparse_array, fnn_create_shortcut_array
  public fnn_free_ann, fnn_randomize_weights, fnn_init_weights
  public fnn_print_connections, fnn_print_parameters, fnn_get_num_input, fnn_get_num_output
  public fnn_get_total_neurons, fnn_get_total_connections, fnn_get_network_type
  public fnn_get_connection_rate, fnn_get_num_layers, fnn_get_layer_array, fnn_get_bias_array
  public fnn_get_connection_array, fnn_set_weight_array, fnn_set_weight
  public fnn_get_decimal_point, fnn_get_multiplier

  public fnn_create_train, fnn_read_train
  public fnn_train, fnn_test, fnn_get_mse, fnn_get_bit_fail, fnn_reset_mse
  public fnn_train_on_data, fnn_train_on_file, fnn_train_epoch, fnn_test_data
  public fnn_read_train_from_file, fnn_create_train_empty, fnn_create_train_from_callback
#ifdef FANN23
  public fnn_create_train_array
  public fnn_get_train_input, fnn_get_train_output, fnn_get_min_train_input
  public fnn_get_max_train_input, fnn_get_min_train_output, fnn_get_max_train_output
#endif
  public fnn_free_train, fnn_shuffle_train_data, fnn_scale_train, fnn_descale_train
  public fnn_set_input_scaling_params, fnn_set_output_scaling_params, fnn_set_scaling_params
  public fnn_clear_scaling_params, fnn_scale_input, fnn_scale_output
  public fnn_descale_input, fnn_descale_output
  public fnn_scale_input_train_data, fnn_scale_output_train_data, fnn_scale_train_data
  public fnn_length_train_data, fnn_num_input_train_data, fnn_num_output_train_data


  public fnn_save_train, fnn_save_train_to_fixed, fnn_save_train_to_file
  public fnn_get_training_algorithm, fnn_set_training_algorithm
  public fnn_get_learning_rate, fnn_set_learning_rate

  public fnn_get_activation_function, fnn_set_activation_function
  public fnn_set_activation_function_layer, fnn_set_activation_function_hidden
  public fnn_set_activation_function_output, fnn_get_activation_steepness
  public fnn_set_activation_steepness, fnn_set_activation_steepness_layer
  public fnn_set_activation_steepness_hidden, fnn_set_activation_steepness_output
  public fnn_get_train_error_function, fnn_set_train_error_function
  public fnn_get_train_stop_function, fnn_set_train_stop_function

  public fnn_get_bit_fail_limit, fnn_set_bit_fail_limit
  public fnn_get_quickprop_decay, fnn_set_quickprop_decay
  public fnn_get_quickprop_mu, fnn_set_quickprop_mu
  public fnn_get_rprop_increase_factor, fnn_set_rprop_increase_factor
  public fnn_get_rprop_decrease_factor, fnn_set_rprop_decrease_factor
  public fnn_get_rprop_delta_max, fnn_set_rprop_delta_max
  public fnn_get_rprop_delta_min, fnn_set_rprop_delta_min
  public fnn_get_rprop_delta_zero, fnn_set_rprop_delta_zero
  public fnn_get_sarprop_weight_decay_shift, fnn_set_sarprop_weight_decay_shift
  public fnn_get_sarprop_step_error_threshold_factor, fnn_set_sarprop_step_error_threshold_factor
  public fnn_get_sarprop_step_error_shift,fnn_set_sarprop_step_error_shift
  public fnn_get_sarprop_temperature, fnn_set_sarprop_temperature

  public fnn_cascade_train
  public fnn_cascadetrain_on_data, fnn_cascadetrain_on_file
  public fnn_get_cascade_output_change_fraction, fnn_set_cascade_output_change_fraction
  public fnn_get_cascade_output_stagnation_epochs, fnn_set_cascade_output_stagnation_epochs
  public fnn_get_cascade_candidate_change_fraction, fnn_set_cascade_candidate_change_fraction
  public fnn_get_cascade_candidate_stagnation_epochs, fnn_set_cascade_candidate_stagnation_epochs
  public fnn_get_cascade_weight_multiplier, fnn_set_cascade_weight_multiplier
  public fnn_get_cascade_candidate_limit, fnn_set_cascade_candidate_limit
  public fnn_get_cascade_max_out_epochs, fnn_set_cascade_max_out_epochs
  public fnn_get_cascade_min_out_epochs, fnn_set_cascade_min_out_epochs
  public fnn_get_cascade_max_cand_epochs, fnn_set_cascade_max_cand_epochs
  public fnn_get_cascade_min_cand_epochs, fnn_set_cascade_min_cand_epochs
  public fnn_get_cascade_num_candidates, fnn_get_cascade_activation_functions_count
  public fnn_get_cascade_activation_functions, fnn_set_cascade_activation_functions
  public fnn_get_cascade_activation_steepnesses_count, fnn_get_cascade_activation_steepnesses
  public fnn_set_cascade_activation_steepnesses, fnn_get_cascade_num_candidate_groups
  public fnn_set_cascade_num_candidate_groups


contains



  subroutine fnn_save(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: ierr
    
    ierr = int(fann_save(ann,f_c_string(filename)))

    if (ierr.ne.0) then
       stop 'fnn_save: error saving neural net!'
    endif

  end subroutine fnn_save



  subroutine fnn_save_to_fixed(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: ierr

    ierr = int(fann_save_to_fixed(ann,f_c_string(filename)))

    if (ierr.ne.0) then
       stop 'fnn_save_to_fixed: error saving neural net!'
    endif

  end subroutine fnn_save_to_fixed




  subroutine fnn_create_from_file(filename)
    implicit none
    character(len=*), intent(in) :: filename

    ann = fann_create_from_file(f_c_string(filename))

  end subroutine fnn_create_from_file




  subroutine fnn_create_standard_array(nlayers,neurons)
    implicit none
    integer(ip), intent(in) :: nlayers
    integer(ip), dimension(nlayers), intent(in) :: neurons

    integer(C_INT) :: num_layers
    integer(C_INT), dimension(nlayers) :: layers

    num_layers = int(nlayers,C_INT)
    layers = int(neurons,C_INT)
       
    ann = fann_create_standard_array(num_layers,layers)
            
  end subroutine fnn_create_standard_array



  subroutine fnn_create_sparse_array(rate,nlayers,neurons)
    implicit none
    real(fp), intent(in) :: rate
    integer(ip), intent(in) :: nlayers
    integer(ip), dimension(nlayers), intent(in) :: neurons

    real(C_FLOAT) :: connection_rate
    integer(C_INT) :: num_layers
    integer(C_INT), dimension(nlayers) :: layers

    connection_rate = real(rate, C_FLOAT)
    num_layers = int(nlayers,C_INT)
    layers = int(neurons,C_INT)
       
    ann = fann_create_sparse_array(connection_rate,num_layers,layers)
            
  end subroutine fnn_create_sparse_array



  subroutine fnn_create_shortcut_array(nlayers,neurons)
    implicit none
    integer(ip), intent(in) :: nlayers
    integer(ip), dimension(nlayers), intent(in) :: neurons

    integer(C_INT) :: num_layers
    integer(C_INT), dimension(nlayers) :: layers

    num_layers = int(nlayers,C_INT)
    layers = int(neurons,C_INT)
       
    ann = fann_create_shortcut_array(num_layers,layers)
            
  end subroutine fnn_create_shortcut_array



  function fnn_run(input)
    implicit none
    real(fp), dimension(:), intent(in) :: input    
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(:), pointer :: fnn_run
    integer(FANN_TYPE), dimension(size(input,1)) :: fannInput
    fannInput = int(input,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(:), pointer :: fnn_run
    real(FANN_TYPE), dimension(size(input,1)) :: fannInput
    fannInput = real(input,FANN_TYPE)
#endif

    fnn_run => f_fann_run(ann,fannInput)        
    
  end function fnn_run


  subroutine fnn_free_ann()
    implicit none
   
    call fann_destroy(ann)

  end subroutine fnn_free_ann

  

  subroutine fnn_print_connections()
    implicit none

    call fann_print_connections(ann)
    
  end subroutine fnn_print_connections


  subroutine fnn_print_parameters()
    implicit none

    call fann_print_parameters(ann)
    
  end subroutine fnn_print_parameters



  function fnn_get_num_input()
    implicit none
    integer(ip) :: fnn_get_num_input

    fnn_get_num_input = int(fann_get_num_input(ann),ip)

  end function fnn_get_num_input



  function fnn_get_num_output()
    implicit none
    integer(ip) :: fnn_get_num_output

    fnn_get_num_output = int(fann_get_num_output(ann),ip)

  end function fnn_get_num_output


     
  function fnn_get_total_neurons()
    implicit none
    integer(ip) :: fnn_get_total_neurons

    fnn_get_total_neurons = int(fann_get_total_neurons(ann),ip)

  end function fnn_get_total_neurons
  


  function fnn_get_total_connections()
    implicit none
    integer(ip) :: fnn_get_total_connections

    fnn_get_total_connections = int(fann_get_total_connections(ann),ip)

  end function fnn_get_total_connections

  

  function fnn_get_network_type()
    implicit none
    integer(ip) :: fnn_get_network_type

    fnn_get_network_type = int(fann_get_network_type(ann),ip)

  end function fnn_get_network_type

  
       
  function fnn_get_connection_rate()
    implicit none
    real(fp) :: fnn_get_connection_rate

    fnn_get_connection_rate = real(fann_get_connection_rate(ann),fp)

  end function fnn_get_connection_rate



  function fnn_get_num_layers()
    implicit none
    integer(ip) :: fnn_get_num_layers

    fnn_get_num_layers = int(fann_get_num_layers(ann),ip)

  end function fnn_get_num_layers
    


  subroutine fnn_get_layer_array(layers)
    implicit none
    integer(ip), dimension(:), intent(inout) :: layers
    integer(C_INT), dimension(size(layers,1)) :: fannLayers
    
    if ((size(layers,1)).ne.fnn_get_num_layers()) then
       stop 'fnn_get_layer_array: layers array size mismatch!'
    endif

    call fann_get_layer_array(ann,fannLayers)

    layers = int(fannLayers,ip)

  end subroutine fnn_get_layer_array



  subroutine fnn_get_bias_array(biases)
    implicit none
    integer(ip), dimension(:), intent(inout) :: biases
    integer(C_INT), dimension(size(biases,1)) :: fannBias

    if ((size(biases,1)).ne.fnn_get_num_layers()) then
       stop 'fnn_get_layer_biases: biases array size mismatch!'
    endif

    call fann_get_bias_array(ann,fannBias)

    biases = int(fannBias,ip)

  end subroutine fnn_get_bias_array



  subroutine fnn_get_connection_array(connections)
    implicit none
    type(fann_connection), dimension(:), intent(inout) :: connections

    if ((size(connections,1)).ne.fnn_get_total_connections()) then
       stop 'fnn_get_connection_array: connections array size mismatch!'
    endif

    call fann_get_connection_array(ann,connections)

  end subroutine fnn_get_connection_array




  subroutine fnn_set_weight_array(connections)
    implicit none
    type(fann_connection), dimension(:), intent(in) :: connections
    integer(C_INT) :: num_connections

    num_connections = size(connections,1)

    if (num_connections.ne.fnn_get_total_connections()) then
       stop 'fnn_set_weight_array: connections array size mismatch!'
    endif

    call fann_set_weight_array(ann,connections,num_connections)

  end subroutine fnn_set_weight_array



  subroutine fnn_set_weight(nfrom, nto, weight)
    implicit none
    integer(ip), intent(in) :: nfrom, nto
    real(fp), intent(in) :: weight
    integer(C_INT) :: from_neuron, to_neuron
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: fannWeight
    fannWeight = int(weight,FANN_TYPE)
#else
    real(FANN_TYPE) :: fannWeight
    fannWeight = real(weight,FANN_TYPE)
#endif
    
    from_neuron = int(nfrom,C_INT)
    to_neuron = int(nto,C_INT)

    call fann_set_weight(ann,from_neuron,to_neuron,fannWeight)   

  end subroutine fnn_set_weight



   function fnn_get_decimal_point()
     implicit none
     integer(ip) :: fnn_get_decimal_point

#ifdef FIXEDFANN     
     fnn_get_decimal_point = int(fann_get_decimal_point(ann),ip)
#else
     stop 'fnn_get_decimal_point: only available for FIXEDFANN!'
#endif

   end function fnn_get_decimal_point



   function fnn_get_multiplier()
     implicit none
     integer(ip) :: fnn_get_multiplier
#ifdef FIXEDFANN     
     fnn_get_multiplier = int(fann_get_multiplier(ann),ip)
#else
     stop 'fnn_get_multiplier: only available for FIXEDFANN!'
#endif

   end function fnn_get_multiplier


  

  subroutine fnn_randomize_weights(min,max)
    implicit none
    real(fp), intent(in) :: min, max
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: min_weight, max_weight
    min_weight = int(min,FANN_TYPE)
    max_weight = int(max,FANN_TYPE)
#else
    real(FANN_TYPE) :: min_weight, max_weight
    min_weight = real(min,FANN_TYPE)
    max_weight = real(max,FANN_TYPE)
#endif
    
    call fann_randomize_weights(ann,min_weight,max_weight)

  end subroutine fnn_randomize_weights


  subroutine fnn_init_weights()
    implicit none

    call fann_init_weights(ann,traindata)

  end subroutine fnn_init_weights



  subroutine fnn_train(intrain,outtrain)
    implicit none
    real(fp), dimension(:), intent(in) :: intrain, outtrain
    logical :: invalid
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(intrain,1)) :: input
    integer(FANN_TYPE), dimension(size(outtrain,1)) :: output
    input = int(intrain,FANN_TYPE)
    output = int(outtrain,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(intrain,1)) :: input
    real(FANN_TYPE), dimension(size(outtrain,1)) :: output
    input = real(intrain,FANN_TYPE)
    output = real(outtrain,FANN_TYPE)
#endif

    invalid = (size(intrain,1).ne.fnn_get_num_input()) &
         .or. (size(outtrain,1).ne.fnn_get_num_output())

    if (invalid) stop 'fnn_train_once: data size does not fit the network!'

    call fann_train(ann,input,output)
    
  end subroutine fnn_train



  function fnn_test(intest,outtest)
    implicit none
    real(fp), dimension(:), intent(in) :: intest
    real(fp), dimension(:), intent(in) :: outtest
    type(C_PTR) :: fannOut
#ifdef FIXEDFANN
    integer(ip), dimension(:), pointer :: fnn_test
    integer(FANN_TYPE), dimension(size(intest,1)) :: fannIn
    integer(FANN_TYPE), dimension(size(outtest,1)) :: fannDesired
    fannIn = int(intest,FANN_TYPE)
    fannDesired = int(outtest,FANN_TYPE)
#else
    real(fp), dimension(:), pointer :: fnn_test
    real(FANN_TYPE), dimension(size(intest,1)) :: fannIn
    real(FANN_TYPE), dimension(size(outtest,1)) :: fannDesired
    fannIn = real(intest,FANN_TYPE)
    fannDesired = real(outtest,FANN_TYPE)
#endif

    fannOut = fann_test(ann,fannIn,fannDesired)

    call C_F_POINTER(fannOut,fnn_test, shape=[fann_get_num_output(ann)])
        
  end function fnn_test



  function fnn_get_mse()
    implicit none
    real(fp) :: fnn_get_mse

    fnn_get_mse = real(fann_get_mse(ann),fp)

  end function fnn_get_mse



  function fnn_get_bit_fail()
    implicit none
    integer(ip) :: fnn_get_bit_fail

    fnn_get_bit_fail = int(fann_get_bit_fail(ann),ip)

  end function fnn_get_bit_fail



  subroutine fnn_reset_mse()
    implicit none

    call fann_reset_mse(ann)

  end subroutine fnn_reset_mse
          


  subroutine fnn_train_on_data(maxEpoch,epochBetweenReport,maxError)
    implicit none
    integer(ip), intent(in) :: maxEpoch, epochBetweenReport
    real(fp), intent(in) :: maxError

    integer(C_INT) :: max_epochs, epochs_between_reports
    real(C_FLOAT) :: desired_error

    max_epochs = int(maxEpoch,C_INT)
    epochs_between_reports = int(epochBetweenReport,C_INT)
    desired_error = real(maxError,C_FLOAT)

    call fann_train_on_data(ann,traindata,max_epochs,epochs_between_reports,desired_error)

  end subroutine fnn_train_on_data




  subroutine fnn_train_on_file(filename,maxEpoch,epochBetweenReport,maxError)
    implicit none
    character(len=*), intent(in) :: filename
    integer(ip), intent(in) :: maxEpoch, epochBetweenReport
    real(fp), intent(in) :: maxError
    
    integer(C_INT) :: max_epochs, epochs_between_reports
    real(C_FLOAT) :: desired_error    

    max_epochs = int(maxEpoch,C_INT)
    epochs_between_reports = int(epochBetweenReport,C_INT)
    desired_error = real(maxError,C_FLOAT)

    call fann_train_on_file(ann,f_c_string(filename),max_epochs,epochs_between_reports,desired_error)

  end subroutine fnn_train_on_file



  function fnn_train_epoch()
    implicit none
    real(fp) :: fnn_train_epoch
    real(C_FLOAT) :: mse

    mse = fann_train_epoch(ann,traindata)

    fnn_train_epoch = real(mse,fp)

  end function fnn_train_epoch



  function fnn_test_data()
    implicit none
    real(fp) :: fnn_test_data
    real(C_FLOAT) :: mse

    mse = fann_test_data(ann,traindata)

    fnn_test_data = real(mse,fp)

  end function fnn_test_data



  subroutine fnn_read_train_from_file(filename)
    implicit none
    character(len=*), intent(in) :: filename

    traindata = fann_read_train_from_file(f_c_string(filename))

  end subroutine fnn_read_train_from_file



  subroutine fnn_create_train_empty(ndata,ndim,nout)
    implicit none
    integer(ip), intent(in) :: ndata, ndim, nout

    integer(C_INT) :: num_data, num_input, num_output

    num_data = int(ndata,C_INT)
    num_input = int(ndim,C_INT)
    num_output = int(nout,C_INT)

    traindata = fann_create_train(num_data, num_input, num_output)

  end subroutine fnn_create_train_empty

  

  subroutine fnn_create_train_from_callback(ndata,ndim,nout,xdata,fdata)
    implicit none
    integer(ip), intent(in) :: ndim, ndata, nout
    real(fp), dimension(ndim,ndata), intent(in) :: xdata
    real(fp), dimension(nout,ndata), intent(in) :: fdata

    integer(C_INT) :: num_data, num_input, num_output

    num_data = int(ndata,C_INT)
    num_input = int(ndim,C_INT)
    num_output = int(nout,C_INT)

    allocate(inTrainData(ndim,ndata))
    allocate(outTrainData(nout,ndata))

    inTrainData(:,:) = xdata(:,:)
    outTrainData(:,:) = fdata(:,:)

    traindata = fann_create_train_from_callback(num_data,num_input,num_output &
         ,C_FUNLOC(fnn_train_callback))

    deallocate(inTrainData,outTrainData)

  end subroutine fnn_create_train_from_callback

  


  subroutine fnn_train_callback(num, num_input, num_output, input, output) bind(C)
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
    
  end subroutine  fnn_train_callback


  
  subroutine fnn_free_train()
    implicit none

    call fann_destroy_train(traindata)

  end subroutine fnn_free_train
  



 subroutine fnn_shuffle_train_data()
    implicit none

    call fann_shuffle_train_data(traindata)

  end subroutine fnn_shuffle_train_data



#ifdef FANN23


  subroutine fnn_create_train_array(ndata,ndim,input,nout,output)
    import C_PTR, FANN_TYPE
    implicit none
    integer(ip), intent(in) :: ndata, ndim, nout
    integer(C_INT) :: num_data, num_input, num_output
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(nin, ndata), intent(in), target :: input
    integer(FANN_TYPE), dimension(nout,ndata), intent(in), target :: output
    integer(FANN_TYPE), dimension(:), pointer :: ptrin => null()
    integer(FANN_TYPE), dimension(:), pointer :: ptrout => null()
#else
    real(FANN_TYPE), dimension(nin,ndata), intent(in), target :: input
    real(FANN_TYPE), dimension(nout,ndata), intent(in), target :: output
    real(FANN_TYPE), dimension(:), pointer :: ptrin => null()
    real(FANN_TYPE), dimension(:), pointer :: ptrout => null()
#endif


    num_data = int(ndata,C_INT)
    num_input = int(ndim,C_INT)
    num_output = int(nout,C_INT)

    !maybe it is the transpose?
    ptrin(0:ndata*ndim-1) => input(1:nin,1:ndata)
    ptrout(0:ndata*nout-1) => output(1:nout,1:ndata)

    traindata = fann_create_train_array(num_data,num_input,ptrin,num_output,ptrout)

    ptrin => null()
    ptrout => null()
  end subroutine fnn_create_train_array


  function fnn_get_train_input(index)
    implicit none
    real(fp) :: fnn_get_train_input
    integer(ip), intent(in) :: index
    integer(C_INT) :: position

    position = int(index,C_INT)

    fnn_get_train_input = real(fann_get_train_input(traindata,position),fp)

  end function fnn_get_train_input



  function fnn_get_train_output(index)
    implicit none
    real(fp) :: fnn_get_train_output
    integer(ip), intent(in) :: index
    integer(C_INT) :: position

    position = int(index,C_INT)

    fnn_get_train_output = real(fann_get_train_output(traindata,position),fp)

  end function fnn_get_train_output


 

  function fnn_get_min_train_input()
    implicit none
    real(fp) :: fnn_get_min_train_input

    fnn_get_min_train_input = real(fann_get_min_train_input(traindata),fp)
    
  end function fnn_get_min_train_input



  function fnn_get_max_train_input()
    implicit none
    real(fp) :: fnn_get_max_train_input

    fnn_get_max_train_input = real(fann_get_max_train_input(traindata),fp)

  end function fnn_get_max_train_input

    

  function fnn_get_min_train_output()
    implicit none
    real(fp) :: fnn_get_min_train_output

    fnn_get_min_train_output = real(fann_get_min_train_output(traindata),fp)
    
  end function fnn_get_min_train_output



  function fnn_get_max_train_output()
    implicit none
    real(fp) :: fnn_get_max_train_output

    fnn_get_max_train_output = real(fann_get_max_train_output(traindata),fp)

  end function fnn_get_max_train_output
#endif


  subroutine fnn_scale_train()
    implicit none

    call fann_scale_train(ann,traindata)

  end subroutine fnn_scale_train




  subroutine fnn_descale_train()
    implicit none

    call fann_descale_train(ann,traindata)

  end subroutine fnn_descale_train


  
  subroutine fnn_set_input_scaling_params(newMin,newMax)
    implicit none
    real(fp), intent(in) :: newMin, newMax
    integer :: ierr
    real(C_FLOAT) :: new_input_min,new_input_max

    new_input_min = real(newMin,C_FLOAT)
    new_input_max = real(newMax,C_FLOAT)

    ierr = int(fann_set_input_scaling_params(ann,traindata,new_input_min,new_input_max))
    
    if (ierr.ne.0) then
       stop 'fnn_set_input_scaling_params: failed!'
    endif

  end subroutine fnn_set_input_scaling_params




  subroutine fnn_set_output_scaling_params(newMin,newMax)
    implicit none
    real(fp), intent(in) :: newMin,newMax
    integer :: ierr
    real(C_FLOAT) :: new_output_min,new_output_max

    new_output_min = real(newMin,C_FLOAT)
    new_output_max = real(newMax,C_FLOAT)

    ierr = int(fann_set_output_scaling_params(ann,traindata,new_output_min,new_output_max))

    if (ierr.ne.0) then
       stop 'fnn_set_output_scaling_params: failed!'
    endif

  end subroutine fnn_set_output_scaling_params




  subroutine fnn_set_scaling_params(newInMin,newInMax,newOutMin,newOutMax)
    implicit none
    real(fp), intent(in) :: newInMin,newInMax,newOutMin,newOutMax
    integer :: ierr
    real(C_FLOAT) :: new_input_min, new_input_max
    real(C_FLOAT) :: new_output_min, new_output_max

    new_input_min = real(newInMin,C_FLOAT)
    new_input_max = real(newInMax,C_FLOAT)
    new_output_min = real(newOutMin,C_FLOAT)
    new_output_max = real(newOutMax,C_FLOAT)

    ierr = int(fann_set_scaling_params(ann,traindata &
         , new_input_min,new_input_max,new_output_min,new_output_max))

    if (ierr.ne.0) then
       stop 'fnn_set_scaling_params: failed!'
    endif

  end subroutine fnn_set_scaling_params



  subroutine fnn_clear_scaling_params()
    implicit none
    integer :: ierr

    ierr = int(fann_clear_scaling_params(ann))

    if (ierr.ne.0) then
       stop 'fnn_clear_scaling_params: failed!'
    endif

  end subroutine fnn_clear_scaling_params




  subroutine fnn_scale_input(input)
    implicit none
    real(fp), dimension(:), intent(inout) :: input
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(input,1)) :: vector
    vector = int(input,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(input,1)) :: vector
#endif
    vector = real(input,FANN_TYPE)

    call fann_scale_input(ann,vector)

    input = real(vector,fp)
   
  end subroutine fnn_scale_input



  
  subroutine fnn_scale_output(output)
    implicit none
    real(fp), dimension(:), intent(inout) :: output
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(output,1)) :: vector
    vector = int(output,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(output,1)) :: vector
#endif
    vector = real(output,FANN_TYPE)

    call fann_scale_output(ann,vector)

    output = real(vector,fp)

  end subroutine fnn_scale_output



  subroutine fnn_descale_input(input)
    implicit none
    real(fp), dimension(:), intent(inout) :: input
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(input,1)) :: vector
    vector = int(input,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(input,1)) :: vector
#endif
    vector = real(input,FANN_TYPE)

    call fann_descale_input(ann,vector)

    input = real(vector,fp)
   
  end subroutine fnn_descale_input



  
  subroutine fnn_descale_output(output)
    implicit none
    real(fp), dimension(:), intent(inout) :: output
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(output,1)) :: vector
    vector = int(output,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(output,1)) :: vector
#endif
    vector = real(output,FANN_TYPE)

    call fann_descale_output(ann,vector)

    output = real(vector,fp)

  end subroutine fnn_descale_output



  subroutine fnn_scale_input_train_data(newMin,newMax)
    implicit none
    real(fp), intent(in) :: newMin, newMax
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: new_min, new_max
    new_min = int(newMin,FANN_TYPE)
    new_max = int(newMax,FANN_TYPE)
#else
    real(FANN_TYPE) :: new_min, new_max
    new_min = real(newMin,FANN_TYPE)
    new_max = real(newMax,FANN_TYPE)
#endif
    
    call fann_scale_input_train_data(traindata,new_min,new_max)

  end subroutine fnn_scale_input_train_data




  subroutine fnn_scale_output_train_data(newMin,newMax)
    implicit none
    real(fp), intent(in) :: newMin, newMax
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: new_min, new_max
    new_min = int(newMin,FANN_TYPE)
    new_max = int(newMax,FANN_TYPE)
#else
    real(FANN_TYPE) :: new_min, new_max
    new_min = real(newMin,FANN_TYPE)
    new_max = real(newMax,FANN_TYPE)
#endif

    call fann_scale_output_train_data(traindata,new_min,new_max)

  end subroutine fnn_scale_output_train_data




  subroutine fnn_scale_train_data(newMin,newMax)
    implicit none
    real(fp), intent(in) :: newMin, newMax
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: new_min, new_max
    new_min = int(newMin,FANN_TYPE)
    new_max = int(newMax,FANN_TYPE)
#else
    real(FANN_TYPE) :: new_min, new_max
    new_min = real(newMin,FANN_TYPE)
    new_max = real(newMax,FANN_TYPE)
#endif

    call fann_scale_train_data(traindata,new_min,new_max)

  end subroutine fnn_scale_train_data




  function fnn_length_train_data()
    implicit none
    integer(ip) :: fnn_length_train_data

    fnn_length_train_data = int(fann_length_train_data(traindata),ip)

  end function fnn_length_train_data
  


  function fnn_num_input_train_data()
    implicit none
    integer(ip) :: fnn_num_input_train_data

    fnn_num_input_train_data = int(fann_num_input_train_data(traindata),ip)

  end function fnn_num_input_train_data
  


  function fnn_num_output_train_data()
    implicit none
    integer(ip) :: fnn_num_output_train_data

    fnn_num_output_train_data = int(fann_num_output_train_data(traindata),ip)

  end function fnn_num_output_train_data


 
  subroutine fnn_save_train_to_file(filename)
    implicit none
    character(len=*), intent(in) :: filename
   
    call fann_save_train(traindata,f_c_string(filename))

  end subroutine fnn_save_train_to_file



  subroutine fnn_save_train_to_fixed(filename,decpoint)
    implicit none
    character(len=*), intent(in) :: filename
    integer(ip), intent(in) :: decpoint    
    integer(C_INT) :: decimal_point

    integer(ip) :: i,n, ierr

    decimal_point = int(decpoint,C_INT)
    
    ierr = int(fann_save_train_to_fixed(traindata,f_c_string(filename),decimal_point))

    if (ierr.ne.0) then
       stop 'fnn_save_train_to_fixed: error in saving file!'
    endif

  end subroutine fnn_save_train_to_fixed



  function fnn_get_training_algorithm()
    implicit none
    integer(ip) :: fnn_get_training_algorithm
    
    fnn_get_training_algorithm = int(fann_get_training_algorithm(ann),ip)

  end function fnn_get_training_algorithm
    


  subroutine fnn_set_training_algorithm(algoname)
    implicit none
    character(len=*), intent(in) :: algoname
    integer(C_INT) :: fannalgo

    fannalgo = enum_training_algorithm(algoname)
       
    call fann_set_training_algorithm(ann,fannalgo)

  end subroutine fnn_set_training_algorithm



  function fnn_get_learning_rate()
    implicit none
    real(fp) :: fnn_get_learning_rate

    fnn_get_learning_rate = real(fann_get_learning_rate(ann),fp)

  end function fnn_get_learning_rate



  subroutine fnn_set_learning_rate(rate)
    implicit none
    real(fp), intent(in) :: rate
    real(C_FLOAT) :: learning_rate

    learning_rate = real(rate,C_FLOAT)

    call fann_set_learning_rate(ann,learning_rate)

  end subroutine fnn_set_learning_rate



  function fnn_get_activation_function(ilayer,ineuron)
    implicit none
    integer(ip) :: fnn_get_activation_function
    integer(ip), intent(in) :: ilayer,ineuron
    integer(C_INT) :: layer, neuron

    layer = int(ilayer,C_INT)
    neuron = int(ineuron,C_INT)

    fnn_get_activation_function = int(fann_get_activation_function(ann,layer,neuron),ip)
    
  end function fnn_get_activation_function



  subroutine fnn_set_activation_function(actiname,ilayer,ineuron)
    implicit none
    character(len=*), intent(in) :: actiname
    integer(ip), intent(in) :: ilayer, ineuron

    integer(C_INT) :: actifunc,layer,neuron

    layer = int(ilayer,C_INT)
    neuron = int(ineuron,C_INT)

    actifunc = enum_activation_function(actiname)

    call fann_set_activation_function(ann,actifunc,layer,neuron)

  end subroutine fnn_set_activation_function



  subroutine fnn_set_activation_function_layer(actiname,ilayer)
    implicit none
    character(len=*), intent(in) :: actiname
    integer(ip), intent(in) :: ilayer
    integer(C_INT) :: layer, actifunc

    layer = int(ilayer,C_INT)

    actifunc = enum_activation_function(actiname)

    call fann_set_activation_function_layer(ann,actifunc,layer)

  end subroutine fnn_set_activation_function_layer



  subroutine fnn_set_activation_function_hidden(actiname)
    implicit none
    character(len=*), intent(in) :: actiname
    integer(C_INT) :: actifunc

    actifunc = enum_activation_function(actiname)

    call fann_set_activation_function_hidden(ann,actifunc)

  end subroutine fnn_set_activation_function_hidden




  subroutine fnn_set_activation_function_output(actiname)
    implicit none
    character(len=*), intent(in) :: actiname
    integer(C_INT) :: actifunc

    actifunc = enum_activation_function(actiname)

    call fann_set_activation_function_output(ann,actifunc)

  end subroutine fnn_set_activation_function_output



  function fnn_get_activation_steepness(ilayer,ineuron)
    implicit none
    real(fp) :: fnn_get_activation_steepness
    integer(ip), intent(in) :: ilayer,ineuron

    integer(C_INT) :: layer, neuron

    layer = int(ilayer,C_INT)
    neuron = int(ineuron,C_INT)

    fnn_get_activation_steepness = real(fann_get_activation_steepness(ann,layer,neuron),fp)

  end function fnn_get_activation_steepness



  subroutine fnn_set_activation_steepness(slope,ilayer,ineuron)
    implicit none
    real(fp), intent(in) :: slope
    integer(ip), intent(in) :: ilayer,ineuron
    integer(C_INT) :: layer, neuron
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: steepness
    steepness = int(slope,FANN_TYPE)
#else
    real(FANN_TYPE) :: steepness
    steepness = real(slope,FANN_TYPE)
#endif

    layer = int(ilayer,C_INT)
    neuron = int(ineuron,C_INT)

    call fann_set_activation_steepness(ann,steepness,layer,neuron)

  end subroutine  fnn_set_activation_steepness



  subroutine fnn_set_activation_steepness_layer(slope,ilayer)
    implicit none
    real(fp), intent(in) :: slope
    integer(ip), intent(in) :: ilayer
    integer(C_INT) :: layer
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: steepness
    steepness = int(slope,FANN_TYPE)
#else
    real(FANN_TYPE) :: steepness
    steepness = real(slope,FANN_TYPE)
#endif

    layer = int(ilayer,C_INT)

    call fann_set_activation_steepness_layer(ann,steepness,layer)

  end subroutine fnn_set_activation_steepness_layer




  subroutine fnn_set_activation_steepness_hidden(slope)
    implicit none
    real(fp), intent(in) :: slope
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: steepness
    steepness = int(slope,FANN_TYPE)
#else
    real(FANN_TYPE) :: steepness
    steepness = real(slope,FANN_TYPE)
#endif

    call fann_set_activation_steepness_hidden(ann,steepness)

  end subroutine fnn_set_activation_steepness_hidden




  subroutine fnn_set_activation_steepness_output(slope)
    implicit none
    real(fp), intent(in) :: slope
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: steepness
    steepness = int(slope,FANN_TYPE)
#else
    real(FANN_TYPE) :: steepness
    steepness = real(slope,FANN_TYPE)
#endif

    call fann_set_activation_steepness_output(ann,steepness)

  end subroutine fnn_set_activation_steepness_output




  function fnn_get_train_error_function()
    implicit none
    integer(ip) :: fnn_get_train_error_function

    fnn_get_train_error_function = int(fann_get_train_error_function(ann),ip)

  end function fnn_get_train_error_function


  subroutine fnn_set_train_error_function(errname)
    implicit none
    character(len=*), intent(in) :: errname
    integer(C_INT) :: errfunc

    errfunc = enum_error_function(errname)

    call fann_set_train_error_function(ann,errfunc)

  end subroutine fnn_set_train_error_function




  function fnn_get_train_stop_function()
    implicit none
    integer(ip) :: fnn_get_train_stop_function

    fnn_get_train_stop_function = int(fann_get_train_stop_function(ann),ip)

  end function fnn_get_train_stop_function


  subroutine fnn_set_train_stop_function(errname)
    implicit none
    character(len=*), intent(in) :: errname
    integer(C_INT) :: errfunc

    errfunc = enum_stop_function(errname)

    call fann_set_train_stop_function(ann,errfunc)

  end subroutine fnn_set_train_stop_function



  function fnn_get_bit_fail_limit()
    implicit none
    real(fp) :: fnn_get_bit_fail_limit

    fnn_get_bit_fail_limit = real(fann_get_bit_fail_limit(ann),fp)

  end function fnn_get_bit_fail_limit



  subroutine fnn_set_bit_fail_limit(limit)
    implicit none
    real(fp), intent(in) :: limit
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: bflimit
    bflimit = int(limit,FANN_TYPE)
#else
    real(FANN_TYPE) :: bflimit
    bflimit = real(limit,FANN_TYPE)
#endif

    call fann_set_bit_fail_limit(ann,bflimit)

  end subroutine fnn_set_bit_fail_limit
    


  function fnn_get_quickprop_decay()
    implicit none
    real(fp) :: fnn_get_quickprop_decay

    fnn_get_quickprop_decay = real(fann_get_quickprop_decay(ann),fp)

  end function fnn_get_quickprop_decay



  subroutine fnn_set_quickprop_decay(decay)
    implicit none
    real(fp), intent(in) :: decay
    real(C_FLOAT) :: quickprop_decay

    quickprop_decay = real(decay,C_FLOAT)

    call fann_set_quickprop_decay(ann,quickprop_decay)

  end subroutine fnn_set_quickprop_decay
  


  function fnn_get_quickprop_mu()
    implicit none
    real(fp) :: fnn_get_quickprop_mu

    fnn_get_quickprop_mu = real(fann_get_quickprop_mu(ann),fp)

  end function fnn_get_quickprop_mu



  subroutine fnn_set_quickprop_mu(mu)
    implicit none
    real(fp), intent(in) :: mu
    real(C_FLOAT) :: quickprop_mu

    quickprop_mu = real(mu,C_FLOAT)

    call fann_set_quickprop_mu(ann,quickprop_mu)

  end subroutine fnn_set_quickprop_mu

    
    
  function fnn_get_rprop_increase_factor()
    implicit none
    real(fp) :: fnn_get_rprop_increase_factor

    fnn_get_rprop_increase_factor = real(fann_get_rprop_increase_factor(ann),fp)
    
  end function fnn_get_rprop_increase_factor


  subroutine fnn_set_rprop_increase_factor(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: rprop_increase_factor

    rprop_increase_factor = real(factor,C_FLOAT)

    call fann_set_rprop_increase_factor(ann,rprop_increase_factor)

  end subroutine fnn_set_rprop_increase_factor


  
  function fnn_get_rprop_decrease_factor()
    implicit none
    real(fp) :: fnn_get_rprop_decrease_factor

    fnn_get_rprop_decrease_factor = real(fann_get_rprop_decrease_factor(ann),fp)

  end function fnn_get_rprop_decrease_factor



  subroutine fnn_set_rprop_decrease_factor(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: rprop_decrease_factor

    rprop_decrease_factor = real(factor,C_FLOAT)
    
    call fann_set_rprop_decrease_factor(ann,rprop_decrease_factor)

  end subroutine fnn_set_rprop_decrease_factor

    

  function fnn_get_rprop_delta_min()
    implicit none
    real(fp) :: fnn_get_rprop_delta_min

    fnn_get_rprop_delta_min = real(fann_get_rprop_delta_min(ann),fp)

  end function fnn_get_rprop_delta_min



  subroutine fnn_set_rprop_delta_min(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: rprop_delta_min

    rprop_delta_min = real(factor,C_FLOAT)

    call fann_set_rprop_delta_min(ann,rprop_delta_min)

  end subroutine fnn_set_rprop_delta_min



  function fnn_get_rprop_delta_max()
    implicit none
    real(fp) :: fnn_get_rprop_delta_max

    fnn_get_rprop_delta_max = real(fann_get_rprop_delta_max(ann),fp)

  end function fnn_get_rprop_delta_max



  subroutine fnn_set_rprop_delta_max(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: rprop_delta_max

    rprop_delta_max = real(factor,C_FLOAT)

    call fann_set_rprop_delta_max(ann,rprop_delta_max)

  end subroutine fnn_set_rprop_delta_max

   
  function fnn_get_rprop_delta_zero()
    implicit none
    real(fp) :: fnn_get_rprop_delta_zero

    fnn_get_rprop_delta_zero = real(fann_get_rprop_delta_zero(ann),fp)

  end function fnn_get_rprop_delta_zero


  subroutine fnn_set_rprop_delta_zero(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: rprop_delta_zero

    rprop_delta_zero = real(factor,C_FLOAT)

    call fann_set_rprop_delta_zero(ann,rprop_delta_zero)

  end subroutine fnn_set_rprop_delta_zero



  function fnn_get_sarprop_weight_decay_shift()
    implicit none
    real(fp) :: fnn_get_sarprop_weight_decay_shift

    fnn_get_sarprop_weight_decay_shift = real(fann_get_sarprop_weight_decay_shift(ann),fp)

  end function fnn_get_sarprop_weight_decay_shift



  subroutine fnn_set_sarprop_weight_decay_shift(shift)
    implicit none
    real(fp), intent(in) :: shift
    real(C_FLOAT) :: sarprop_weight_decay_shift

    sarprop_weight_decay_shift = real(shift,C_FLOAT)

    call fann_set_sarprop_weight_decay_shift(ann,sarprop_weight_decay_shift)

  end subroutine fnn_set_sarprop_weight_decay_shift



  function fnn_get_sarprop_step_error_threshold_factor()
    implicit none
    real(fp) :: fnn_get_sarprop_step_error_threshold_factor

    fnn_get_sarprop_step_error_threshold_factor &
         = real(fann_get_sarprop_step_error_threshold_factor(ann),fp)

  end function fnn_get_sarprop_step_error_threshold_factor


  subroutine fnn_set_sarprop_step_error_threshold_factor(factor)
    implicit none
    real(fp), intent(in) :: factor
    real(C_FLOAT) :: sarprop_step_error_threshold_factor

    sarprop_step_error_threshold_factor = real(factor,C_FLOAT)

    call fann_set_sarprop_step_error_threshold_factor(ann,sarprop_step_error_threshold_factor)

  end subroutine fnn_set_sarprop_step_error_threshold_factor



    
  function fnn_get_sarprop_step_error_shift()
    implicit none
    real(fp) :: fnn_get_sarprop_step_error_shift

    fnn_get_sarprop_step_error_shift = real(fann_get_sarprop_step_error_shift(ann),fp)

  end function fnn_get_sarprop_step_error_shift
  
  
  

  subroutine fnn_set_sarprop_step_error_shift(shift)
    implicit none
    real(fp), intent(in) :: shift
    real(C_FLOAT) :: sarprop_step_error_shift

    sarprop_step_error_shift = real(shift,C_FLOAT)

    call fann_set_sarprop_step_error_shift(ann,sarprop_step_error_shift)

  end subroutine fnn_set_sarprop_step_error_shift




  function fnn_get_sarprop_temperature()
    implicit none
    real(fp) :: fnn_get_sarprop_temperature

    fnn_get_sarprop_temperature = real(fann_get_sarprop_temperature(ann),fp)

  end function fnn_get_sarprop_temperature



  subroutine fnn_set_sarprop_temperature(temp)
    implicit none
    real(fp), intent(in) :: temp
    real(C_FLOAT) :: sarprop_temperature

    sarprop_temperature = real(temp,C_FLOAT)

    call fann_set_sarprop_temperature(ann,sarprop_temperature)

  end subroutine fnn_set_sarprop_temperature



  subroutine fnn_cascadetrain_on_data(maxNeurons,neuronsBetweenReports,maxError)
    implicit none
    integer(ip), intent(in) :: maxNeurons, neuronsBetweenReports
    real(fp), intent(in) :: maxError

    integer(C_INT) :: max_neurons, neurons_between_reports
    real(C_FLOAT) :: desired_error

    max_neurons = int(maxNeurons,C_INT)
    neurons_between_reports = int(neuronsBetweenReports,C_INT)
    desired_error = real(maxError,C_FLOAT)
  
    call fann_cascadetrain_on_data(ann,traindata, max_neurons, neurons_between_reports, desired_error)
  
  end subroutine fnn_cascadetrain_on_data



  subroutine fnn_cascadetrain_on_file(filename,maxNeurons,neuronsBetweenReports,maxError)
    implicit none
    character(len=*), intent(in) :: filename
    integer(ip), intent(in) :: maxNeurons, neuronsBetweenReports
    real(fp), intent(in) :: maxError

    integer(C_INT) :: max_neurons, neurons_between_reports
    real(C_FLOAT) :: desired_error

    max_neurons = int(maxNeurons,C_INT)
    neurons_between_reports = int(neuronsBetweenReports,C_INT)
    desired_error = real(maxError,C_FLOAT)
  
    call fann_cascadetrain_on_file(ann, f_c_string(filename), max_neurons, neurons_between_reports, desired_error)
  
  end subroutine fnn_cascadetrain_on_file




  function fnn_get_cascade_output_change_fraction()
    implicit none
    real(fp) :: fnn_get_cascade_output_change_fraction

    fnn_get_cascade_output_change_fraction = real(fann_get_cascade_output_change_fraction(ann),fp)

  end function fnn_get_cascade_output_change_fraction



  subroutine fnn_set_cascade_output_change_fraction(fraction)
    implicit none
    real(fp), intent(in) :: fraction
    real(C_FLOAT) :: change_fraction

    change_fraction = real(fraction, C_FLOAT)

    call fann_set_cascade_output_change_fraction(ann,change_fraction)

  end subroutine fnn_set_cascade_output_change_fraction



  function fnn_get_cascade_output_stagnation_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_output_stagnation_epochs

    fnn_get_cascade_output_stagnation_epochs = int(fann_get_cascade_output_stagnation_epochs(ann),ip)

  end function fnn_get_cascade_output_stagnation_epochs



  subroutine fnn_set_cascade_output_stagnation_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: stagnation_epochs

    stagnation_epochs = int(epochs, C_INT)

    call fann_set_cascade_output_stagnation_epochs(ann,stagnation_epochs)

  end subroutine fnn_set_cascade_output_stagnation_epochs



  function fnn_get_cascade_candidate_change_fraction()
    implicit none
    real(fp) :: fnn_get_cascade_candidate_change_fraction

    fnn_get_cascade_candidate_change_fraction = real(fann_get_cascade_candidate_change_fraction(ann),fp)

  end function fnn_get_cascade_candidate_change_fraction



  subroutine fnn_set_cascade_candidate_change_fraction(fraction)
    implicit none
    real(fp), intent(in) :: fraction
    real(C_FLOAT) :: change_fraction

    change_fraction = real(fraction, C_FLOAT)

    call fann_set_cascade_candidate_change_fraction(ann,change_fraction)

  end subroutine fnn_set_cascade_candidate_change_fraction




  function fnn_get_cascade_candidate_stagnation_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_candidate_stagnation_epochs

    fnn_get_cascade_candidate_stagnation_epochs = int(fann_get_cascade_candidate_stagnation_epochs(ann),ip)

  end function fnn_get_cascade_candidate_stagnation_epochs



  subroutine fnn_set_cascade_candidate_stagnation_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: stagnation_epochs

    stagnation_epochs = int(epochs, C_INT)

    call fann_set_cascade_candidate_stagnation_epochs(ann,stagnation_epochs)

  end subroutine fnn_set_cascade_candidate_stagnation_epochs

    


  function fnn_get_cascade_weight_multiplier()
    implicit none
    real(fp) :: fnn_get_cascade_weight_multiplier

    fnn_get_cascade_weight_multiplier = real(fann_get_cascade_weight_multiplier(ann),fp)

  end function fnn_get_cascade_weight_multiplier




  subroutine fnn_set_cascade_weight_multiplier(multiplier)
    implicit none
    real(fp), intent(in) :: multiplier
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: weight_multiplier
    weight_multiplier = int(multiplier,FANN_TYPE)
#else
    real(FANN_TYPE) :: weight_multiplier
    weight_multiplier = real(multiplier,FANN_TYPE)
#endif
    call fann_set_cascade_weight_multiplier(ann,weight_multiplier)

  end subroutine fnn_set_cascade_weight_multiplier




  function fnn_get_cascade_candidate_limit()
    implicit none
    real(fp) :: fnn_get_cascade_candidate_limit

    fnn_get_cascade_candidate_limit = real(fann_get_cascade_candidate_limit(ann),fp)

  end function fnn_get_cascade_candidate_limit




  subroutine fnn_set_cascade_candidate_limit(limit)
    implicit none
    real(fp), intent(in) :: limit
#ifdef FIXEDFANN
    integer(FANN_TYPE) :: candidate_limit
    candidate_limit = int(limit,FANN_TYPE)
#else
    real(FANN_TYPE) :: candidate_limit
    candidate_limit = real(limit,FANN_TYPE)
#endif
    call fann_set_cascade_candidate_limit(ann,candidate_limit)

  end subroutine fnn_set_cascade_candidate_limit

    

  
  function fnn_get_cascade_max_out_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_max_out_epochs

    fnn_get_cascade_max_out_epochs = int(fann_get_cascade_max_out_epochs(ann),ip)

  end function fnn_get_cascade_max_out_epochs



  subroutine fnn_set_cascade_max_out_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: max_out_epochs

    max_out_epochs = int(epochs,C_INT)

    call fann_set_cascade_max_out_epochs(ann,max_out_epochs)

  end subroutine fnn_set_cascade_max_out_epochs
  

  function fnn_get_cascade_min_out_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_min_out_epochs

    fnn_get_cascade_min_out_epochs = int(fann_get_cascade_min_out_epochs(ann),ip)

  end function fnn_get_cascade_min_out_epochs



  subroutine fnn_set_cascade_min_out_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: min_out_epochs

    min_out_epochs = int(epochs,C_INT)

    call fann_set_cascade_min_out_epochs(ann,min_out_epochs)

  end subroutine fnn_set_cascade_min_out_epochs



  function fnn_get_cascade_max_cand_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_max_cand_epochs

    fnn_get_cascade_max_cand_epochs = int(fann_get_cascade_max_cand_epochs(ann),ip)

  end function fnn_get_cascade_max_cand_epochs



  subroutine fnn_set_cascade_max_cand_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: max_cand_epochs

    max_cand_epochs = int(epochs,C_INT)

    call fann_set_cascade_max_cand_epochs(ann,max_cand_epochs)

  end subroutine fnn_set_cascade_max_cand_epochs



  function fnn_get_cascade_min_cand_epochs()
    implicit none
    integer(ip) :: fnn_get_cascade_min_cand_epochs

    fnn_get_cascade_min_cand_epochs = int(fann_get_cascade_min_cand_epochs(ann),ip)

  end function fnn_get_cascade_min_cand_epochs



  subroutine fnn_set_cascade_min_cand_epochs(epochs)
    implicit none
    integer(ip), intent(in) :: epochs
    integer(C_INT) :: min_cand_epochs

    min_cand_epochs = int(epochs,C_INT)

    call fann_set_cascade_min_cand_epochs(ann,min_cand_epochs)

  end subroutine fnn_set_cascade_min_cand_epochs


     

  function fnn_get_cascade_num_candidates()
    implicit none
    integer(ip) :: fnn_get_cascade_num_candidates

    fnn_get_cascade_num_candidates = int(fann_get_cascade_num_candidates(ann),ip)

  end function fnn_get_cascade_num_candidates




  function fnn_get_cascade_activation_functions_count()
    implicit none
    integer(ip) :: fnn_get_cascade_activation_functions_count

    fnn_get_cascade_activation_functions_count = int(fann_get_cascade_activation_functions_count(ann),ip)

  end function fnn_get_cascade_activation_functions_count


    

  function fnn_get_cascade_activation_functions()
    implicit none
    integer(C_INT), dimension(:), pointer :: fnn_get_cascade_activation_functions
    
    fnn_get_cascade_activation_functions => f_fann_get_cascade_activation_functions(ann)
    
  end function fnn_get_cascade_activation_functions




  subroutine fnn_set_cascade_activation_functions(actinames)
    implicit none
    character(len=*), dimension(:), intent(in) :: actinames

    integer(C_INT), dimension(size(actinames,1)) :: activation_functions_id
    
    integer :: n,i
    integer(C_INT) :: activation_functions_count

    n = size(actinames,1)
    activation_functions_count = int(n,C_INT)

    do i=1,n
       activation_functions_id(i) = enum_activation_function(trim(actinames(i)))
    end do

    call fann_set_cascade_activation_functions(ann,activation_functions_id,activation_functions_count)
        
  end subroutine fnn_set_cascade_activation_functions




  function fnn_get_cascade_activation_steepnesses_count()
    implicit none
    integer(ip) :: fnn_get_cascade_activation_steepnesses_count
  
    fnn_get_cascade_activation_steepnesses_count &
         = int(fann_get_cascade_activation_steepnesses_count(ann),ip)
  
  end function fnn_get_cascade_activation_steepnesses_count



  function fnn_get_cascade_activation_steepnesses()
    implicit none
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(:), pointer :: fnn_get_cascade_activation_steepnesses
#else
    real(FANN_TYPE), dimension(:), pointer :: fnn_get_cascade_activation_steepnesses
#endif

    fnn_get_cascade_activation_steepnesses => f_fann_get_cascade_activation_steepnesses(ann)

  end function fnn_get_cascade_activation_steepnesses



  subroutine fnn_set_cascade_activation_steepnesses(steeps)
    implicit none
    real(fp), dimension(:), intent(in) :: steeps
    integer(C_INT) :: activation_steepnesses_count
    integer :: n
#ifdef FIXEDFANN
    integer(FANN_TYPE), dimension(size(steeps,1)) :: activation_steepnesses
    activation_steepnesses = int(steeps,FANN_TYPE)
#else
    real(FANN_TYPE), dimension(size(steeps,1)) :: activation_steepnesses
    activation_steepnesses = real(steeps,FANN_TYPE)
#endif

    n = size(steeps,1)

    activation_steepnesses_count = int(n,C_INT)
    
    call fann_set_cascade_activation_steepnesses(ann,activation_steepnesses &
         ,activation_steepnesses_count)
   
  end subroutine fnn_set_cascade_activation_steepnesses




  function fnn_get_cascade_num_candidate_groups()
    implicit none
    integer(ip) :: fnn_get_cascade_num_candidate_groups

    fnn_get_cascade_num_candidate_groups = int(fann_get_cascade_num_candidate_groups(ann),ip)

  end function fnn_get_cascade_num_candidate_groups



  subroutine fnn_set_cascade_num_candidate_groups(ngroups)
    implicit none
    integer(ip), intent(in) :: ngroups
    integer(C_INT) :: num_candidate_groups
    
    num_candidate_groups = int(ngroups,C_INT)

    call fann_set_cascade_num_candidate_groups(ann,num_candidate_groups)
    
  end subroutine fnn_set_cascade_num_candidate_groups

     



end module fnn

