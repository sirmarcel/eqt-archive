# set to some unique directory (useful to choose SLURM job id on HPC to avoid problems with concurrent jobs)
socketdir=/tmp/ipi-$SLURM_JOB_ID

echo "starting..."
i-pi -S $socketdir input-mol.xml | tee ipi.out &
echo "started i-pi"
sleep 2
echo "starting driver..."
# adapt quadrature_order and inversions to particular experiment (this is 2i)
i-pi-py_driver -S $socketdir -u -m pet -a driver -o "../model,init-mol.xyz,device=cuda,checkpoint=model,quadrature_order=2,inversions=True" | tee calc.out &
echo "started everything... let's go"
wait
