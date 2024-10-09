# set to some unique directory (useful to choose SLURM job id on HPC to avoid problems with concurrent jobs)
socketdir=/tmp/ipi-$SLURM_JOB_ID

echo "starting..."
i-pi -S $socketdir input-mts.xml | tee ipi.out &
echo "started i-pi"
sleep 2
echo "starting driver..."
i-pi-py_driver -S $socketdir -u -m pet -a driver -o "$original/../model,init-512_h2o.xyz,device=cuda,checkpoint=model" | tee calc.out &
# adapt quadrature_order and inversions to particular experiment; this is 2i
i-pi-py_driver -S $socketdir -u -m pet -a driver-aug -o "$original/../model,init-512_h2o.xyz,device=cuda,checkpoint=model,quadrature_order=2,inversions=True" | tee calc-aug.out &

echo "started everything... let's go"
wait
