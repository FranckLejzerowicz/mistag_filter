import random
import itertools

pairs={}
for i in open('design.txt'):
    pairs[tuple([i.split(',')[1], i.split(',')[2].strip()])]=1

o=open('fasta_noMistag.fasta', 'w')

for i in ['AAAAAAAAAAAAAAAAAAAAAAAAAA', 'CCCCCCCCCCCCCCCCCCCCCCCC']:
  for s in pairs:
    o.write('>%s;size=%s;fwd=%s;rev=%s\n%s\n' % (i, random.choice(range(1000,100000)), s[0], s[1], i))

f=[]
r=[]
for i in open('design.txt'):
  if 'sample' not in i:
    f.append(i.split(',')[1])
    r.append(i.strip().split(',')[2])

for i in ['AAAAAAAAAAAAAAAAAAAAAAAAAA', 'CCCCCCCCCCCCCCCCCCCCCCCC']:
  for j in itertools.product(list(set(f)),list(set(r))):
    if j in pairs:
      continue
    o.write('>%s;size=%s;fwd=%s;rev=%s\n%s\n' % (i, random.choice(range(1,100)), j[0], j[1], i))

o.close()
