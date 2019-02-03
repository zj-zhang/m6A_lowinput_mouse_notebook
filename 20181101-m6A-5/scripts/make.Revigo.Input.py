import sys
n=0
for l in open(sys.argv[1]):
	if n==0:
		n+=1
		continue
	ls=l.strip().split('\t')
	if float(ls[4])<=0.01:
		print ls[0].strip('"').split(' ')[-1].strip('(').strip(')'), ls[4]
