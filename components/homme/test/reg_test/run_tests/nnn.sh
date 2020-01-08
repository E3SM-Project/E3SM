
source mmm2.sh

sed -e s/BBNAME/${name}/ -e s/BBTTYPE/${ttype}/  \
     theta-xx-test-template.cmake > ${name}.cmake



