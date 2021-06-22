#!/bin/bash

# mail -s "Script $0 started on jasmin" nele.reyniers@live.be

python3 /home/users/nelerey/phd_code/ETo/eto/methods/ETo_chess.py \
-v huss,psl,rls,rss,sfcWind,tas \
-i /home/users/nelerey/data/toydata_3x3/{var}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc \
-o /home/users/nelerey/data/toydata_3x3/pet_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc
