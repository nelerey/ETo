
# SLURM_ARRAY_TASK_ID_sf=`printf %02d $SLURM_ARRAY_TASK_ID`
for start in {1981..2051..10};
do
    stop="$((start + 29))"
#    echo pr_rcp85_land-rcm_uk_12km_${SLURM_ARRAY_TASK_ID_sf}_day_${start}0101-${stop}1230_BC.nc
done

python3 /home/users/nelerey/phd_code/ETo/applications/apply_etofao-56_ukcp18_jasmin-ceda.py \
-v hurs,psl,rls,rss,sfcWind,tasmax,tasmin,tas \
-i /badc/ukcp18/data/land-rcm/uk/12km/rcp85/01/{var}/day/latest/{var}_rcp85_land-rcm_uk_12km_01_day_19801201-19901130.nc \
-o /home/users/nelerey/data/PE/ETo-fao56/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc \
-n evpot-fao56
