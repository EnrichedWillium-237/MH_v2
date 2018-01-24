#!/bin/sh

root -l -q -b 'GetVN.C+("../MH.root","N1MCm22SUB2",0.0,2.4,false,-1,1.2499,-0.05,0.05,-0.1,0.03)'
root -l -q -b 'GetVN.C+("../MH.root","N1MCp22SUB2",-2.4,0.0,false,-1,1.2499,-0.05,0.05,-0.1,0.03)'

root -l -q -b 'GetVN.C+("../MH.root","N1MCm22SUB3",0.0,2.4,false,-1,1.2499,-0.05,0.05,-0.1,0.03)'
root -l -q -b 'GetVN.C+("../MH.root","N1MCp22SUB3",-2.4,0.0,false,-1,1.2499,-0.05,0.05,-0.1,0.03)'
