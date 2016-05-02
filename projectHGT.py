#!/usr/bin/python2.7

'''
things to do:
1) make plot of abundance per Hamm Dist class
2) ancestor's trace
3) multiple deaths + large pop size
4) make multiple runs with variable HGT length
5) ...one day make length itself evolvable
6) make command line save file... and while you are at it also parameters
'''

import sys,getopt, math,os,random,string,math
import operator
import functools

#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

#mutations
POINT_MUT_EVENT=True 
HGT_EVENT=True

#mutations of mutation rates
MUT_COMPETENCE_EVENT=False
MUT_MUTRATE=False

#genome
alphabet=['1','0']	#mutations and initialisation make use of this, this is extendable
lengen=20	#genome length
maxfitness=float(lengen)	#maximum fitness

#initialize global lists
lpop=[]   #population list
toplot=[] #plotting data 
lTime=[]  #also some plotting data

################################
### Read commandline options ###
################################

#default values:
popsize=500
MaxGenerations=5000
maxlenHGT=lengen/2 - 1
fileout="hgt_data.txt"
init_competence=0.25
init_mutrate=0.05

try:
  opts, args = getopt.getopt(sys.argv[1:],"hi:g:l:f:c:m:",["help","popsize=","MaxGen=","maxlenHGT=","fileout=","init_comp=","init_mut="])
except getopt.GetoptError:
  print 'projectHGT.py [--popsize,-i] [--MaxGen,-g] [--maxlenHGT, -l] [--fileout, -f] [--init_comp, -c] [--init_mut, -m]'
  sys.exit(2)

for opt, arg in opts:
  if opt in ("-h", "--help"):
    print 'projectHGT.py [--popsize,-i] [--MaxGen,-g] [--maxlenHGT, -l] [--fileout, -f] [--init_comp, -c] [--init_mut, -m]'
    print 'default values:\n'
    print 'popsize: '+str(popsize)
    print 'MaxGenerations: '+str(MaxGenerations)
    print 'MaxlengthHGT: '+str(maxlenHGT)
    print 'fileout: '+fileout
    print 'init_competence: '+str(init_competence)
    print 'init_mutationrate: '+str(init_mutrate)

    print '\nMutations:'
    print 'point mutations: '+str(POINT_MUT_EVENT)
    print 'HGT mutations: '+str(HGT_EVENT)
    print 'Mut. of competence: '+str(MUT_COMPETENCE_EVENT)
    print 'Mut. of mutrate: '+str(MUT_MUTRATE)
    sys.exit()
  elif opt in ("-i", "--popsize"):
    popsize=int(arg)
  elif opt in ("--MaxGen","-g"):
    MaxGenerations=int(arg)
  elif opt in ("--maxlenHGT", "-l"):
    maxlenHGT=int(arg)
  elif opt in ("--fileout", "-f"):
    fileout=arg
  elif opt in ("--init_comp", "-c"):
    init_competence=float(arg)
  elif opt in ("--init_mut", "-m"):
    init_mutrate=float(arg)
  else:
    print 'Unrecognised option, ignored.'


#derived parameters
MaxTime=MaxGenerations*popsize
if maxlenHGT==0: HGT_EVENT=False

############################
### Save options to file ###
############################

#now that we do not have the actual values in a file, it is handy to store the used options in a file for later reference
parfile=fileout+'.par'
pardata=open(parfile,'w')
pardata.write("point mutations: "+str(POINT_MUT_EVENT)+"\n")
pardata.write("HGT mutations: "+str(HGT_EVENT)+"\n")
pardata.write("mutation of mutrate: "+str(MUT_MUTRATE)+"\n")
pardata.write("mutation of HGT competence: "+str(MUT_COMPETENCE_EVENT)+"\n")
pardata.write("\nGenome:\n")
pardata.write("alphabet used: "+str(alphabet)+"\n")
pardata.write("genome length: "+str(lengen)+"\n")
pardata.write("maxfitness: "+str(maxfitness)+"\n")
pardata.write("\ncommandline:\n")

for opt, arg in opts:
  pardata.write(opt+" "+arg+"\n")
pardata.close()

#################
### Functions ###
#################

def Fitness(genome):
  score=sum( [int(pos) for pos in genome ] )
  #score=functools.reduce(operator.mul,[1.5*int(pos) for pos in genome ] , 1)
  #print genome,score
  return score/maxfitness

#takes a genome and randomly mutates the position rpos 
# into a different base randomly chosen from alphabet
def PointMut(genome):
  rpos=random.randrange(lengen)
  bla=list(alphabet)
  bla.remove( genome[rpos] )
  newbase=random.choice(bla)
  newgenome= genome[:rpos] + newbase + genome[rpos+1:]
  
  return newgenome
  
############################
##### --- PRINTING --- #####
############################

def PrintEntropy(Time, lpop):
  lS=(lengen+1)*[0]
  for ic in lpop:
    pos=int(float(lengen)*ic[1])
    lS[pos]+=1
    
  fpopsize=float(popsize)
  
  #in theory there are 2^(lengen) different genomes, 
  #in practice there can be at most 500, i.e. pop_size
  genomes={}
  for ic in lpop:
    if ic[0] in genomes: genomes[ic[0]]+=1
    else: genomes[ic[0]]=1
  
  #print genomes
  entropy=-sum( [ ( genomes[x] /fpopsize)*math.log(genomes[x]/fpopsize) for x in genomes ] )
  
  '''
  bestfit=-1.
  for ic in lpop:
   if ic[1]>bestfit: 
     bestfitgenome=ic[0]
     bestfit=ic[1]
     bestcomp=ic[2]
  '''
  
  avcomp=sum( x[2] for x in lpop )/fpopsize
  avmutrate=sum( x[3] for x in lpop )/fpopsize
  #print "Time",Time,"Entropy",entropy,"Fittest guy",bestfit
  print "Time",Time, "D",lS, "av cmp", avcomp, "av mut",avmutrate, "S", entropy
    
  toplot.append([Time,[x/fpopsize for x in lS], avcomp,avmutrate,entropy])
  #print toplot, lS
  #for x,y in zip(toplot,lS):
  #  print x,y
  #  x.append(y)
    
    
  #print Time, toplot
  
##################################
##### --- INITIALISATION --- #####
##################################

for i in range(popsize):
  #genome= ''.join(random.choice(bases) for _ in range(lengen))
  genome=''.join('0' for _ in range(lengen))
  genome=PointMut(genome)	# we make 1 random mutation per genome
  fitness_genome=Fitness(genome)
  lpop.append([ genome , fitness_genome ,init_competence,init_mutrate ])
  #print Fitness(lpop[i][0])

#############################
##### --- MAIN LOOP --- #####
#############################
  
for Time in range(MaxTime):
  
  if Time%1000==0:
    PrintEntropy(Time, lpop)
  
  
  tot_fitness=sum( [ic[1] for ic in lpop] )
  random_fitness=random.uniform(0, tot_fitness)
  current=0.
  for i,ic in enumerate(lpop):
    current += ic[1]
    if current > random_fitness: 
      #this is the guy that reproduces
      
      #this is where it reproduces
      rand_child_pos=random.randrange(popsize)
      lpop[rand_child_pos]=list(ic)
      
      #MUTATIONS - point mutations
      if POINT_MUT_EVENT==True and random.random()< lpop[rand_child_pos][3]:
        lpop[rand_child_pos][0]=PointMut(lpop[rand_child_pos][0]) #Point mutations
        lpop[rand_child_pos][1]=Fitness(lpop[rand_child_pos][0]) #Fitness of the mutant
        
      #MUTATIONS - competence
      if MUT_COMPETENCE_EVENT==True and random.random()< lpop[rand_child_pos][3]: 
        lpop[rand_child_pos][2] += random.uniform(-0.025, 0.025)
        if lpop[rand_child_pos][2]<0.: lpop[rand_child_pos][2] = -lpop[rand_child_pos][2]
        if lpop[rand_child_pos][2]>1.: lpop[rand_child_pos][2] = 2.-lpop[rand_child_pos][2]
        
      #MUTATIONS - mutation rate
      if MUT_MUTRATE==True and random.random()< lpop[rand_child_pos][3]: 
        lpop[rand_child_pos][3] += random.uniform(-0.025, 0.025)
        if lpop[rand_child_pos][3]<0.: lpop[rand_child_pos][3] = -lpop[rand_child_pos][3]
        if lpop[rand_child_pos][3]>1.: lpop[rand_child_pos][3] = 2.-lpop[rand_child_pos][3]
        
      #MUTATIONS - HGT
      if HGT_EVENT==True and random.random()< lpop[rand_child_pos][2]:
        ric=random.choice(lpop) #pick random individual
        if ric[0]==lpop[rand_child_pos][0]: 
          break
        
        lenHGT=random.randrange(maxlenHGT) #random length
        rpos=random.randrange(lengen) #random position for insertion
        if rpos+lenHGT>=lengen: lenHGT=lengen-rpos
        if ric[0][rpos:rpos+lenHGT] == lpop[rand_child_pos][0][rpos:rpos+lenHGT]: break
        lpop[rand_child_pos][0]= lpop[rand_child_pos][0][:rpos]+ric[0][rpos:rpos+lenHGT] +lpop[rand_child_pos][0][rpos+lenHGT:]
        
      break;

#the list of stuff to print
#Time, lS, avcomp, avmutrate, entropy

toplot=zip(*toplot)
lTime=toplot[0]

#print len(toplot),len(lTime),len(toplot[1])

nsubplots=3
# Two subplots, unpack the axes array immediately
fig, ax = plt.subplots(nsubplots, 1, sharex=True)

#prints frequencies
subplotnum=0

cm = plt.get_cmap('Dark2')
tpl=zip(*toplot[1])
ax[subplotnum].set_color_cycle([cm(1.*i/(lengen+1)) for i in range(lengen+1)])

for i in range(lengen+1):
    ax[subplotnum].plot(lTime, tpl[i], label=str(20-i) )

#for i,ll in enumerate( zip(*toplot[1]) ):
#  ax[subplotnum].plot(lTime,ll,label=str(20-i))
  
  
subplotnum=1
ax[subplotnum].plot(lTime, toplot[2] ,label="av compet")
ax[subplotnum].plot(lTime, toplot[3] ,label="av mutrate")
  
subplotnum=2
ax[subplotnum].plot(lTime, toplot[4] ,label="entropy")

fontsize="x-small"
for n in range(nsubplots):
  # Shrink current axis by 20%
  box = ax[n].get_position()
  ax[n].set_position([box.x0, box.y0, box.width * 0.85, box.height])
  # Put a legend to the right of the current axis
  ax[n].legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fontsize,ncol=3)

plt.savefig("HGT_in_time.png")
toplot=zip(*toplot)
with open(fileout,'w') as fout:
  for line in toplot:
    bla=[fout.write(str(x)+' ') for x in line]
    fout.write('\n')
    
#plt.show()  














