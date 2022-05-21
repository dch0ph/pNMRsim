	function cosh (exp(#1)+exp(-#1))*0.5
	function sinh (exp(#1)-exp(-#1))*0.5 
	function tanh_ (#1-(1/#1))/(#1+(1/#1))
	function tanh tanh_(exp(#1))
	function log10 ln(#1)/ln(10)
	function norm if(size(#1),norm(tail(#1))+head(#1)^2,0) L
