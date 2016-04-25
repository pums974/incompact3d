#
set terminal jpeg size 1080,720 # Sortie jpeg
set output ’tampon.jpeg’ # Nom du fichier de sortie
set nokey # Sans titre
set pm3d map # Trace d’une carte
set palette rgbformulae 33,13,10 # Degrade bleu-rouge
set cbrange[-1:1] # Etendue de la coloration
set tics out # Tics vers l’exterieur
unset ztics # Pas de graduation selon z
splot ’tampon.dat’
#
set output ’temp.jpeg’ # Nom du fichier de sortie
set nokey # Sans titre
set pm3d map # Trace d’une carte
set palette rgbformulae 33,13,10 # Degrade bleu-rouge
set cbrange[292:293] # Etendue de la coloration
set tics out # Tics vers l’exterieur
unset ztics # Pas de graduation selon z
splot ’temp.dat’
