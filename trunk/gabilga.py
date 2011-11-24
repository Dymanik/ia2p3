from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import G1DBinaryString
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Mutators
from pyevolve import DBAdapters
import csv

import sys

from random import randint as rand_randint, choice as rand_choice
from random import random as rand_random
from random import shuffle

from pyevolve import Util

#probability to apply add alternative
PADD=0.01
#probability to apply drop condition
PDROP=0.60

ruleSize=41
atribNum=10
# v1 -16 17-20 21-24 25-28 29-32 33-36 37-40 41-44 45-48 49+
# v4 0 1 2 3 4 5-8 9+
atribSize=[10,4 ,4 ,7 ,2 ,2 ,4 ,4 ,2 ,2 ]
atribPos= [0 ,10,14,18,25,27,29,33,37,39]

trainset=[]
inputset=[]
valset=[]

	
def evaluate(genome, x):
	xt=[]
	if x[0]<=16 : xt+=[0]
	elif 17<=x[0]<=20: xt+=[1]
	elif 21<=x[0]<=24: xt+=[2]
	elif 25<=x[0]<=28: xt+=[3]
	elif 29<=x[0]<=32: xt+=[4]
	elif 33<=x[0]<=36: xt+=[5]
	elif 37<=x[0]<=40: xt+=[6]
	elif 41<=x[0]<=44: xt+=[7]
	elif 45<=x[0]<=48: xt+=[8]
	else: xt+=[9]

	xt+=[x[1]-1]
	xt+=[x[2]-1]

	if x[3]<=0 : xt+=[0]
	elif x[3]==1: xt+=[1]
	elif x[3]==2: xt+=[2]
	elif x[3]==3: xt+=[3]
	elif x[3]==4: xt+=[4]
	elif 5<=x[3]<=8: xt+=[5]
	else: xt+=[6]

	xt+=[x[4]]
	xt+=[x[5]]
	xt+=[x[6]-1]
	xt+=[x[7]-1]
	xt+=[x[8]]
	xt+=[x[9]]

	rules = len(genome)/ruleSize

	y=[]
	
	for i in xrange(rules):
		value =0
		for j in xrange(ruleSize-2):
			if genome[i*ruleSize+j]==1:value-=1
		for j in xrange(9):
		#	print i*ruleSize,atribPos[j],xt[j]
			if genome[i*ruleSize+atribPos[j]+xt[j]]==1:value+=atribSize[j]
		ans = genome[(i+1)*ruleSize-1]
		ans += genome[(i+1)*ruleSize-2]

		y+=[(value,ans)]

	y.sort(key=lambda a:a[0])
	return y[-1][1]


# This function is the evaluation function, we want
# to give high score to more zero'ed chromosomes
def eval_func(chromosome):
	score = 0.0
#	print trainset
	if len(chromosome) > ruleSize*25:return score
	for x in trainset:
#		print evaluate(chromosome,x), x[9]
		if evaluate(chromosome,x)==x[9]-1:score +=1
	score= score/len(trainset)
	return score**2


def DropCondition(genome, **args):
#	printg(genome)

	stringLength = len(genome)
	mutations = PDROP * (stringLength/ruleSize*(atribNum))
	if mutations < 1.0:
		mutations = 0
		for it in xrange(stringLength/ruleSize*atribNum):
			if Util.randomFlipCoin(PDROP):
				for i in xrange(atribSize[it%atribNum]-1):
					genome[ruleSize*(it/atribNum) + atribPos[it%atribNum] + i]=1
				mutations +=1
	else:
		for it in xrange(int(round(mutations))):
			which = rand_randint(0, stringLength/ruleSize*atribNum-1)
			for i in xrange(atribSize[which%atribNum]):
					genome[ruleSize*(which/atribNum) + atribPos[which%atribNum] + i]=1
#	print '>>', mutations
#	printg (genome)
	return int(mutations)
					



def AddAlternative(genome, **args):
   """ The classical flip mutator for binary strings """
   stringLength = len(genome)
   mutations = PADD * (stringLength/ruleSize*atribNum)
   
   if mutations < 1.0:
      mutations = 0
      for it in xrange(stringLength/ruleSize*atribNum):
         if Util.randomFlipCoin(PADD):
            while True:
                zeros=0
                for i in xrange(atribSize[it%atribNum]):
                    if genome[ruleSize*(it/atribNum) + atribPos[it%atribNum] + i] ==0:
                        zeros+=1
                        if Util.randomFlipCoin(0.5):
                           genome[ruleSize*(it/atribNum) + atribPos[it%atribNum] + i] = 1
                           mutations+=1
                           break
                if mutations >0 or zeros==0:break
   else:
      for it in xrange(int(round(mutations))):
         while True:    
            zeros=0
            m=0
            for i in xrange(atribSize[it%atribNum]):
                if genome[ruleSize*(it/atribNum) + atribPos[it%atribNum] + i] ==0:
                    zeros+=1
                    if Util.randomFlipCoin(0.5):
                       genome[ruleSize*(it/atribNum) + atribPos[it%atribNum] + i] = 1
                       mutations+=1
                       m=1
                       break
            if m >0 or zeros==0:break
   return int(mutations)

def CrossOver(genome,**args):
   """ The 1D Binary String crossover, Two Point

   .. warning:: You can't use this crossover method for binary strings with length of 1.

   """
   sister = None
   brother = None
   gMom = args["mom"]
   gDad = args["dad"]
   
   if len(gMom) == 1:
      Util.raiseException("The Binary String have one element, can't use the Two Point Crossover method !", TypeError)

   if len(gMom)<=len(gDad):
      Mom=gMom
      Dad=gDad
   else:
      Mom=gDad
      Dad=gMom

   cuts = [rand_randint(1, len(Mom)-1), rand_randint(1, len(Mom)-1)]
   
   if cuts[0] > cuts[1]:
      Util.listSwapElement(cuts, 0, 1)

   posible1=[]
   posible2=[]

   for it in xrange(len(Dad)/ruleSize):
      posible1+=[(cuts[0]%ruleSize)+ruleSize*it]
      posible2+=[(cuts[1]%ruleSize)+ruleSize*it]

#   print len(Mom),len(Dad),cuts,posible1,posible2
   posible1 = filter(lambda x:x<=posible2[-1],posible1)
   cuts2=[posible1[rand_randint(0,len(posible1)-1)]]
   t=filter(lambda x:x>=cuts2[0],posible2)

   #print cuts2[0],t
   cuts2+=[t[rand_randint(0,len(t)-1)]]
   #print len(Dad),len(Mom),cuts,posible1,posible2,t,cuts2
   if args["count"] >= 1:
      sister = Mom.clone()
      sister.resetStats()
      sister[cuts[0]:cuts[1]] = Dad[cuts2[0]:cuts2[1]]
      sister.genomeSize = len(sister)

   if args["count"] == 2:
      brother = Dad.clone()
      brother.resetStats()
      brother[cuts2[0]:cuts2[1]] = Mom[cuts[0]:cuts[1]]
      brother.genomeSize = len(brother)

#   print Mom.getBinary(),len(Mom)
#   print Dad.getBinary(),len(Dad)
#   print '>>',cuts,cuts2
#   print sister.getBinary(),len(sister)
#   print brother.getBinary(),len(brother)
   return (sister, brother)



def printg(genome):
	rules = len(genome)/ruleSize
	for x in xrange(rules):
		for i in xrange(10):
			for j in xrange(atribSize[i]):
				sys.stdout.write(str(genome[x*ruleSize+atribPos[i]+j]))
			sys.stdout.write(" ")
		sys.stdout.write(",")
	sys.stdout.write("\n")

def run_main():

   Reader = csv.reader(open('Dataextra.txt','rb'))

   i=0
   for x in Reader:
      inputset.insert(i,map(int,x))
      i+=1

   shuffle(inputset)
   p=0.5
   n= int(round(len(inputset)*p))
 
   trainset[0:n] = inputset[0:n]
#   print trainset
   valset = inputset[n+1:-1]

   # Genome instance
   genome = G1DBinaryString.G1DBinaryString(ruleSize*10)

   # The evaluator function (objective function)
   genome.evaluator.set(eval_func)
   genome.mutator.set(Mutators.G1DBinaryStringMutatorFlip)
   genome.mutator.add(AddAlternative)
   genome.mutator.add(DropCondition)
   genome.mutator.setRandomApply(True)
   genome.crossover.set(CrossOver)

   # Genetic Algorithm Instance
   ga = GSimpleGA.GSimpleGA(genome)
   ga.selector.set(Selectors.GTournamentSelector)
   sqlite_adapter = DBAdapters.DBSQLite(identify=sys.argv[1])
   ga.setDBAdapter(sqlite_adapter)
   ga.setGenerations(70)
   ga.setElitism(True)

   # Do the evolution, with stats dump
   # frequency of 10 generations
   ga.evolve(freq_stats=0)

   # Best individual
   #print ga.bestIndividual()
   #printg (ga.bestIndividual())

   score = 0.0
   for x in valset:
      if evaluate(ga.bestIndividual(),x)==x[9]-1:score +=1

   score = score/len(valset)
   print score
	

if __name__ == "__main__":
   run_main()
   
