from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import G1DBinaryString
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Mutators
import csv

from random import randint as rand_randint, choice as rand_choice
from random import random as rand_random
from random import shuffle

from pyevolve import Util

ruleSize=76
atribNum=10
atribSize=[34,4,4,17,2,2,4,4,2,3]

trainset=[]
inputset=[]
valset=[]

# This function is the evaluation function, we want
# to give high score to more zero'ed chromosomes
def eval_func(chromosome):
   score = 0.0

   for i in xrange(len(chromosome)):
       if chromosome.genomeList[i]==1:score+=1
      
   return score/len(chromosome)

def AddAlternative(genome, **args):
   """ The classical flip mutator for binary strings """
   if args["pmut"] <= 0.0: return 0
   stringLength = len(genome)
   mutations = args["pmut"] * (stringLength/ruleSize*atribNum)
   
   if mutations < 1.0:
      mutations = 0
      for it in xrange(stringLength/ruleSize*atribNum):
         if Util.randomFlipCoin(args["pmut"]):
            while True:
                zeros=0
                for i in xrange(atribSize[it%atribNum]):
                    if genome[it/atribNum + it%atribNum +i ] ==0:
                        zeros+=1
                        if Util.randomFlipCoin(0.5):
                           genome[it/atribNum + it%atribNum +i ] = 1
                           mutations+=1
                           break
                if mutations >0 or zeros==0:break
   else:
      for it in xrange(int(round(mutations))):
         while True:    
            zeros=0
            m=0
            for i in xrange(atribSize[it%atribNum]):
                if genome[it/atribNum + it%atribNum +i ] ==0:
                    zeros+=1
                    if Util.randomFlipCoin(0.5):
                       genome[it/atribNum + it%atribNum +i ] = 1
                       mutations+=1
                       m=1
                       break
            if m >0 or zeros==0:break
   return int(mutations)

def AddAlternative(genome, **args):
   """ The classical flip mutator for binary strings """
   if args["pmut"] <= 0.0: return 0
   stringLength = len(genome)
   mutations = args["pmut"] * (stringLength/ruleSize*atribNum)
   
   if mutations < 1.0:
      mutations = 0
      for it in xrange(stringLength/ruleSize*atribNum):
         if Util.randomFlipCoin(args["pmut"]):
            while True:
                zeros=0
                for i in xrange(atribSize[it%atribNum]):
                    if genome[it/atribNum + it%atribNum +i ] ==0:
                        zeros+=1
                        if Util.randomFlipCoin(0.5):
                           genome[it/atribNum + it%atribNum +i ] = 1
                           mutations+=1
                           break
                if mutations >0 or zeros==0:break
   else:
      for it in xrange(int(round(mutations))):
         while True:    
            zeros=0
            m=0
            for i in xrange(atribSize[it%atribNum]):
                if genome[it/atribNum + it%atribNum +i ] ==0:
                    zeros+=1
                    if Util.randomFlipCoin(0.5):
                       genome[it/atribNum + it%atribNum +i ] = 1
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



def run_main():

   Reader = csv.reader(open('Dataextra.txt','rb'))

   i=0
   for x in Reader:
      inputset.insert(i,x)
      i+=1

   shuffle(inputset)
   p=0.5
   n= int(round(len(inputset)*p))   
   trainset = inputset[0:n]
   valset = inputset[n+1:-1]

   # Genome instance
   genome = G1DBinaryString.G1DBinaryString(ruleSize*2)

   # The evaluator function (objective function)
   genome.evaluator.set(eval_func)
   genome.mutator.set(Mutators.G1DBinaryStringMutatorFlip)
   genome.crossover.set(CrossOver)

   # Genetic Algorithm Instance
   ga = GSimpleGA.GSimpleGA(genome)
   ga.selector.set(Selectors.GRouletteWheel)
   ga.setGenerations(70)

   # Do the evolution, with stats dump
   # frequency of 10 generations
   ga.evolve(freq_stats=20)

   # Best individual
   print ga.bestIndividual()

if __name__ == "__main__":
   run_main()
   
