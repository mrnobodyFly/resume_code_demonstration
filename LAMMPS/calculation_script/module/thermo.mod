thermo_style custom	step temp epair pe ke etotal press vol &
pxx pyy pzz pxy pxz pyz lx ly lz xlo xhi ylo yhi zlo zhi fnorm fmax
thermo_modify	lost error flush yes format 4 %.8f
thermo		${nthermo}
