#!/bin/bash

# A. Bullones
# 2/07/22
# last modification 16/11/22
# merge all the files obtained with FeatureCounts

# eliminar el archivo mergefile.txt por si estuviera creado
# para que no duplique información
rm counts_m.txt

# para quedarnos con las columnas que nos interesan
for num in *.txt;
do cut -f 1,6 $num > $num.cut;
done

# para renombrar las columnas
for num in *.cut;
do
  name=${num/\.txt\.cut/}
  # eliminamos la cabecera
  tail -n +3 $num > $num.middle
  # añadimos la nueva cabecera
  echo -e "ID\t$name" | cat - $num.middle > $num.name

  rm $num.middle
done

# eliminar los archivos .cut
rm *.cut

# unir todos los archivos en uno
awk 'NF > 0 { a[$1] = a[$1] "\t" $2 } END { for (i in a) { print i a[i]; } }' *cut.name > length_m.txt

# ponemos la cabecera como primera línea
# guardamos la cabecera en una variable
header=$(grep 'ID' length_m.txt)
# eliminamos la cabecera
sed -i '/'"$header"'/d' length_m.txt
# la añadimos al principio
sed -i '1i '"$header"'' length_m.txt

# eliminar los archivos cut.name
rm *cut.name

echo "done :)"
