optopt { 
	$1 $(3?$2)
	autoopt 1e-4 $(3?:$1)
}
