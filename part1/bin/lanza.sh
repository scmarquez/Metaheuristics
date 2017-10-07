#!/bin/bash
archivos=(chr20b.dat chr22a.dat els19.dat esc32b.dat kra30b.dat lipa90b.dat nug25.dat sko56.dat sko64.dat sko72.dat sko100a.dat sko100b.dat sko100c.dat sko100d.dat sko100e.dat tai30b.dat tai50b.dat tai60a.dat tai256c.dat tho150.dat)

for (( i = 0; i < 20; i++ )); do
	#statements
	./p1 ${archivos[i]} 3 >> salida.txt
done

