#!/bin/bash

for f in frame-*; do
  cp $f init.xyz
  i=${f#frame-}
  i=${i%.xyz}
  if [ ! -e ice-geop-${i}_aug.out ]; then  
    echo Running $i
    i-pi input-geop.xml &
    sleep 2
#    i-pi-py_driver -u -m pet -a driver -o "../model,init.xyz,device=cpu,checkpoint=model,quadrature_order=2,inversions=True" &> log.pet-aug
    i-pi-py_driver -u -m pet -a driver -o "../model,init.xyz,device=cpu,checkpoint=model" &> log.pet
    wait
    for a in ice-geop.*; do
        mv $a ice-geop-${i}.${a#ice-geop.}
    done
    echo Running $i with augmentation
    i-pi input-geop.xml &
    sleep 2
    i-pi-py_driver -u -m pet -a driver -o "../model,init.xyz,device=cpu,checkpoint=model,quadrature_order=2,inversions=True" &> log.pet-aug
    wait
    for a in ice-geop.*; do
        mv $a ice-geop-${i}_aug.${a#ice-geop.}
    done
  fi
done
    
