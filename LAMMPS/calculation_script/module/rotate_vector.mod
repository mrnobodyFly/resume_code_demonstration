# rotate vector n1 n2 n3
# ouput n1g n2g n3g
variable n1g equal "v_R11 * v_n1 + v_R12 * v_n2 + v_R13 * v_n3"
variable n2g equal "v_R21 * v_n1 + v_R22 * v_n2 + v_R23 * v_n3" 
variable n3g equal "v_R31 * v_n1 + v_R32 * v_n2 + v_R33 * v_n3" 

if "${n1g} > 0" then &
   "variable n1g equal -${n1g}" &
   "variable n2g equal -${n2g}" &
   "variable n3g equal -${n3g}"
