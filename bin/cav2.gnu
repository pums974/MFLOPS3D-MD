
file0="fort.20"
file1="fort.21"
file2="fort.22"
file3="fort.23"
file4="fort.24"
file5="fort.25"
file6="fort.26"
file7="fort.27"
file8="fort.28"
file9="fort.29"



set title "Erreur vitesse u"
splot file0 u 1:2:3 w lp,file1 u 1:2:3 w lp,file2 u 1:2:3 w lp,file3 u 1:2:3 w lp
pause -1

set title "Erreur vitesse v"
splot file0 u 1:2:4 w lp,file1 u 1:2:4 w lp,file2 u 1:2:4 w lp,file3 u 1:2:4 w lp
pause -1

set title "Erreur vitesse w"
splot file0 u 1:2:5 w lp,file1 u 1:2:5 w lp,file2 u 1:2:5 w lp,file3 u 1:2:5 w lp
pause -1

set title "Erreur vitesse p"
splot file0 u 1:2:6 w lp,file1 u 1:2:6 w lp,file2 u 1:2:6 w lp,file3 u 1:2:6 w lp
pause -1

set title "Erreur vitesse div"
splot file0 u 1:2:7 w lp,file1 u 1:2:7 w lp,file2 u 1:2:7 w lp,file3 u 1:2:7 w lp
pause -1

