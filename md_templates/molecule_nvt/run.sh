# set to some unique directory (useful to choose SLURM job id on HPC to avoid problems with concurrent jobs)
socketdir=/tmp/ipi-$SLURM_JOB_ID

echo "starting..."
i-pi -S $socketdir input-mol.xml | tee ipi.out &
echo "started i-pi"
sleep 2
echo "starting driver..."
i-pi-py_driver -S $socketdir -u -m pet -a driver -o "../model,init-mol.xyz,device=cuda,checkpoint=model" | tee calc-og.out &
i-pi-py_driver -S $socketdir -u -m pet -a driver-aug -o "../model,init-mol.xyz,device=cuda,checkpoint=model,quadrature_order=2,inversions=True" | tee calc-aug.out &
echo "started everything... let's go"
wait
