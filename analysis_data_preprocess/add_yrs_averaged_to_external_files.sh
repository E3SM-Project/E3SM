path='/p/cscratch/acme/data/obs_for_acme_diags_v1.6.0'
cd $path

declare -a sn=("ANN" "DJF" "MAM" "JJA" "SON")
for j in "${sn[@]}"
do
    ncatted -O -a yrs_averaged,global,a,c,"198307-200806" ISCCPCOSP*${j}_climo.nc
    ncatted -O -a yrs_averaged,global,a,c,"200701-201001" CALIPSOCOSP*${j}_climo.nc
done




declare -a sn=("ANN" "DJF" "JJA")
for j in "${sn[@]}"
do
    ncatted -O -a yrs_averaged,global,a,c,"1987-2000" SSMI_${j}_climo.nc
done
