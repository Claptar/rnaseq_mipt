pdf("RSeQCresults/delition_profile/SRR3414636.deletion_profile.pdf")
pos=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64)
value=c(0,0,0,0,0,0,0,2,429,432,383,571,621,622,871,921,1011,1066,1164,1101,1122,1249,1150,1385,1180,1171,1261,1245,1114,1326,1176,1134,1037,1108,1093,1268,1256,1227,1043,1199,1155,1284,1228,1163,1112,1142,976,981,967,1066,1102,981,999,848,815,631,539,485,2,1,0,0,0,0,0)
plot(pos,value,type='b', col='blue',xlab="Read position (5'->3')", ylab='Deletion count')
dev.off()
