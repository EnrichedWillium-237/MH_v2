#!/bin/sh
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-0.8,0.8,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-2.4,2.4,true)'

root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-2.4,-2.0,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-2.0,-1.6,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-1.2,-0.8,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",-0.8,-0.4,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",0.0,0.4,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",0.4,0.8,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",0.8,1.2,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",1.2,1.6,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",1.6,2.0,true)'
root -l -q -b 'GetVN.C+("../MH.root","N1EVENSUB3",2.0,2.4,true)'


