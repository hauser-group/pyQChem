import sys
a=open("ch3sh_nh3_cs_x8.xyz",'r').read().splitlines()
l=12
b=len(a)/l

c=[90,95,100,105,110,125,150,200]

for i in xrange(b):
	o=open("ch3sh_nh3_cs_x8.xyz".replace("x8",str(c[i])),'w')
	for j in xrange(l):
		o.write(a[i*l+j]+'\n')
	
	o.close()
