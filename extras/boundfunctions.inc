# Use 'apply' functions to convert 'external' parameter to 'internal' parameter that can vary freely
# 'remove' functions reverse the transformation

# P_int = sqrt( (P_ext-a+1)^2-1)  lower bound on P_ext (#1) is a (#2)
	function apply_lower_bound sqrt((#1-#2+1)^2-1)
# P_ext = a-1+sqrt(P_int^2+1)
	function remove_lower_bound #2-1+sqrt(#1^2+1)

# P_int = sqrt( (b-P_ext+1)^2-1)  upper bound on P_ext (#1) is b (#2)
	function apply_upper_bound sqrt((#2-#1+1)^2-1)
# P_ext = b+1-sqrt(P_int^2+1)
	function remove_upper_bound #2+1-sqrt(#2^2+1)

# double-sided bound a (#2) to b (#3)
# N.B. modified from Minuit so that scale is preserved in linear region (around mid-point of a and b)
	function apply_double_bound_ asin((#1-#2)/#3)
	function apply_double_bound apply_double_bound_(#1,0.5*(#2+#3),0.5*(#3-#2))
	function remove_double_bound_ #2+#3*sin(#1)
	function remove_double_bound remove_double_bound_(#1,0.5*(#2+#3),0.5*(#3-#2))
			  