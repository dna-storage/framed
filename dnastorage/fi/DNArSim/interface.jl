#Author: Kevin Volkel
#Description: part of a DNArSim port for FrameD infrastructure


k=0
root=""
#simple hook to bring in python arguments
function load_parameters(arg_k,arg_root)
	 global root="$(arg_root)/k$(arg_k)"
	 global k=arg_k
end