
source create-theta-temporary.sh

sed -e s/BBNAME/${name}/ \
    -e s/BBIFKOKKOS/${ifkokkos}/ \
    -e s/BBADVFORM/${advform}/  \
    -e s/BBTTYPE/${ttype}/  \
    -e s/BBHVS/${hvs}/  \
    -e s/BBTOM/${hvst}/  \
    -e s/BBRSPLIT/${rsplit}/  \
    -e s/BBQSIZE/${qsize}/  \
    -e s/BBNTOP/${nutop}/  \
    -e s/BBNDIV/${nudiv}/  \
     create-theta-template.cmake > ${name}.cmake



