	function factorial if(#1,#1*factorial(#1-1),1)
	function sum if(size(#1),head(#1)+sum(tail(#1)),0) L
	function complement [head(#1,#2-1),tail(#1,#2)] LS
	function rawfourier if(#1,extract(#2,#1)*sin(#1*#3)+rawfourier(#1-1,#2,#3),0) SLL
	function fourier rawfourier(size(#1),#1,#2) LL
	variable result fourier([10 20],[0:10:180])
