import sys
f=open(sys.argv[1],'r')
lines=f.readlines()
f.close()
print(len(lines))
