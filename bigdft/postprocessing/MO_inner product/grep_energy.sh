
#!/bin/bash

#virtuals-k001-NR.b000080
#wavefunction-k001-NR.b000090

name_template=wavefunction-k001-NR.b

num_start=1
num_occ=121
num_vir=10
num_total=208 #$(($num_occ+$num_vir))


for i in $(seq $num_start $num_total);
do
	grep  '# '$(printf %05d $i) log-2ASP-109_E00.yaml>>occupency_energy.log
done

