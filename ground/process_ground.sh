#! /bin/bash

#run GRaWACpy processing.
## specify sys arguments: 
## dates:
yyyy='yyyy'
mm='mm'
dd='dd'
## water vapor retrieval specifications: R: vertical resolution in m, tavg: temporal resolution in s;
R=200
tavg=120

## run processing:
for dd in {23..28}
do
    echo "processing $yyyy $mm $dd"
    python process_grawacpy_ground.py $yyyy $mm $dd $R $tavg
done
