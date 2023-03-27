import os,sys
fh=open(sys.argv[1],'r')

count=0
for line in fh:
 s=line.strip().split('\t')
 count=count+1
 outline1='\t'.join([s[0],"0",s[1],"PEDV_"+str(count),"NCExon","+"])
 outline2='\t'.join([s[0],s[1],s[4],"PEDV_"+str(count),"NCIntron","+"])
 outline3='\t'.join([s[0],s[4],"28044","PEDV_"+str(count),"NCExon","+"])
 print(outline1)
 print(outline2)
 print(outline3)

fh.close()
