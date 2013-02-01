
file0="fort.10"
file1="fort.11"
file2="fort.12"
file3="fort.13"
file4="fort.14"
file5="fort.15"
file6="fort.16"
file7="fort.17"
file8="fort.18"
file9="fort.19"
file10="fort.20"
file11="fort.21"
file12="fort.22"
file13="fort.23"
file14="fort.24"
file15="fort.25"


# 8 procs 10-12-14-16
#file0="fort.10"
#file1="fort.12"
#file2="fort.14"
#file3="fort.16"
#set zrange[-0.1:0.1]

splot file0 u 1:2:3 w lp,file1 u 1:2:3 w lp,file2 u 1:2:3 w lp,file3 u 1:2:3 w lp
pause -1
splot file0 u 1:2:3 w lp,file1 u 1:2:3 w lp,file2 u 1:2:3 w lp,file3 u 1:2:3 w lp,file0 u 1:2:9 w lp,file1 u 1:2:9 w lp,file2 u 1:2:9 w lp,file3 u 1:2:9 w lp
pause -1

splot file0 u 1:2:4 w lp,file1 u 1:2:4 w lp,file2 u 1:2:4 w lp,file3 u 1:2:4 w lp
pause -1
splot file0 u 1:2:4 w lp,file1 u 1:2:4 w lp,file2 u 1:2:4 w lp,file3 u 1:2:4 w lp,file0 u 1:2:10 w lp,file1 u 1:2:10 w lp,file2 u 1:2:10 w lp,file3 u 1:2:10 w lp
pause -1

splot file0 u 1:2:5 w lp,file1 u 1:2:5 w lp,file2 u 1:2:5 w lp,file3 u 1:2:5 w lp
pause -1
splot file0 u 1:2:5 w lp,file1 u 1:2:5 w lp,file2 u 1:2:5 w lp,file3 u 1:2:5 w lp,file0 u 1:2:11 w lp,file1 u 1:2:11 w lp,file2 u 1:2:11 w lp,file3 u 1:2:11 w lp
pause -1

#set zrange[-0.01:0.06]
splot file0 u 1:2:6 w lp,file1 u 1:2:6 w lp,file2 u 1:2:6 w lp,file3 u 1:2:6 w lp
pause -1
splot file0 u 1:2:6 w lp,file1 u 1:2:6 w lp,file2 u 1:2:6 w lp,file3 u 1:2:6 w lp,file0 u 1:2:8 w lp,file1 u 1:2:8 w lp,file2 u 1:2:8 w lp,file3 u 1:2:8 w lp
pause -1


splot file0 u 1:2:7 w lp,file1 u 1:2:7 w lp,file2 u 1:2:7 w lp,file3 u 1:2:7 w lp
pause -1

splot file0 u 1:2:7 w lp,file1 u 1:2:7 w lp,file2 u 1:2:7 w lp,file3 u 1:2:7 w lp,file0 u 1:2:8 w lp,file1 u 1:2:8 w lp,file2 u 1:2:8 w lp,file3 u 1:2:8 w lp
pause -1
