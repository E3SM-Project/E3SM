n=$1
rest=($@)
unset rest[0]
rest="${rest[@]}"

if [ $n -lt 6 ]; then
    r=$n
else
    r=6
fi

cmd="jsrun -n $n -r $r -l gpu-gpu -b packed:1 -d plane:1 -a 1 -c 7 -g 1 --smpiargs \"-gpu\" $rest"
echo $cmd
eval $cmd
