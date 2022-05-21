#exact modelling of phase transients

function transx #2+#4*(1-#2+#1*$lambdaQ*#3) LSSL
function transy #3+#4*(-#3+#1*$lambdaQ*(1-#2)) LSSL
function transphase_ atan(transy(#1,#2,#3,#4),transx(#1,#2,#3,#4)) LSSL
function transphase transphase_(#1,cos(#2),sin(#2),exp(-$lambdaT*#1)) LS
function transamp_ sqrt(transx(#1,#2,#3,#4)^2+transy(#1,#2,#3,#4)^2) LSSL
function transamp transamp_(#1,cos(#2),sin(#2),exp(-$lambdaT*#1)) LS
function makescale (#1/#2)*[1:#2]
function expandtransphase_ if(size(#1),[#2+transphase(#3,head(#1)-#2),expandtransphase_(tail(#1),head(#1),#3)],[]) LSL
function expandtransamp_ if(size(#1),[transamp(#3,head(#1)-#2),expandtransamp_(tail(#1),head(#1),#3)],[]) LSL
function expandtransphase expandtransphase_(#1,#2,makescale(#3,#4)) LSSS
function expandtransamp expandtransamp_(#1,#2,makescale(#3,#4)) LSSS
function expandtransphase expandtransphase_(#1,last(#1),makescale(#2,#3)) LSS
function expandtransamp expandtransamp_(#1,last(#1),makescale(#2,#3)) LSS
 
function cumul if(size(#1),[#2+0.5*head(#1),cumul(tail(#1),#2+head(#1))],[]) LS
function expandtrans__ [repeat(#4,size(#2)),#1*expandtransamp_(#2,#3,#5),expandtransphase_(#2,#3,#5)] SLSLL
function expandtrans_ expandtrans__(#1,#2,#3,#4,cumul(#4,0)) SLSL
function maketscale_ repeat(#1/#2,#2)
function maketscale maketscale_(#1,ceil(#1/$microT))
function makeextscale_ if(#2>#1,[maketscale(#1),WarnOnce('Transient response is not complete within pulse duration')],[maketscale(#2),#1-#2])
function makeextscale if($lambdaT,makeextscale_(#1,$timeconstantsT/$lambdaT),#1)
function expandtrans expandtrans_(#1,#2,#3,makeextscale(#4)) SLSS
function expandtrans expandtrans_(#1,#2,last(#2),makeextscale(#3)) SLS

