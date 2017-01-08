#/usr/bin/python
import math
import re
from string import digits
import random
import numpy as np

class makerFunctions:
	def __init__(self,fileName):
		self._fN=fileName
		self._density=3
		self._angles=False
		self._angleString=[]
		self._angleSequence=[]
		self._angleType=[]
		self._numAngles=0
		

		
		self.importData()
		self.systemName()
		self.parseSequence()
		self.calcSystemParameters()
		
		#self._nAngleTypes=len(self._angleType)
		self._nAngleTypes=len(set(self._angleType))
		
		self.printSystemInfo()
		
		#Create arrays for chain coordinates, bond info
		coordData=np.zeros(shape=(int(self._numAtoms),6))
		
	def writeDataFile(self):
		self.initializeDataStructures()
		self.assignCoordinatesAndBonds()
		self._fo=open(self._nameString, "w")
		self.writeHeader()
		self._fo.write("\nAtoms\n")
		self.writeStructuretoFile(self._atomInformation)
		self._fo.write("\n\nBonds\n")
		self.writeStructuretoFile(self._bondInformation)
		self._fo.write("\n\nAngles\n")
		self.writeStructuretoFile(self._angleInformation)
		print self._nameString,"file written\n"
		self._fo.close()

		
	def writeHeader(self):
		self._fo.write("LAMMPS data file for system %s\n\n" %self._nameString)
		##System info
		self._fo.write("%i atoms\n" %self._numAtoms)
		self._fo.write("%s bonds\n" %int(self._numBonds))
		self._fo.write("%s angles\n" %self._numAngles)
		self._fo.write("0 dihedrals\n")
		self._fo.write("0 impropers\n")
		self._fo.write("%s atom types\n" %int(self._numTypes))
		self._fo.write("%s bond types\n" %self._nBondTypes)
		self._fo.write("%s angle types\n\n" %self._nAngleTypes)
		for dim in ["x","y","z"]:
			strDim=dim+"lo"
			strDim2=dim+"hi"
			self._fo.write("0 %0.3f %s %s\n" %(self._boxLength, strDim, strDim2))
			
	def writeStructuretoFile(self, data):
		for item in data:
			if np.any(item)==1:
				self._fo.write("\n")
				for item2 in item:
					num=float(item2)
					if (num.is_integer()):
						num=int(num)    
					self._fo.write("%s " %num)
		
	
	def assignCoordinatesAndBonds(self):
		for ch in range(0,int(self._numChains)):
			startIndex=(ch)*(self._numMonMol)+1
			startBondIndex=ch*self._numBondsPolymer+1
			startAngleIndex=ch*(self._numMonMol-2)+1
			self.chainCoordinates(ch,startIndex)
			self.chainBonds(ch,startBondIndex)
			if (self._angles==True):
				self.chainAngles(ch,startAngleIndex)
		watStart=int(self._numMonMol)*int(self._numChains)+1
		random.seed()
		for n in range(0,int(self._numWat)):
			x=watStart+n
			self._atomInformation[x-1,0]=x
			self._atomInformation[x-1,1]=int(self._numChains)+n+1
			self._atomInformation[x-1,2]=self._numTypes
			self._atomInformation[x-1,3]=random.random()*self._boxLength
			self._atomInformation[x-1,4]=random.random()*self._boxLength
			self._atomInformation[x-1,5]=random.random()*self._boxLength
		
	def importData(self):
		#f=open(self._fN)
		f=self._fN.split("\n")
		for line in f:
			p=line.split('=')
			if (p[0]=='numtypes'):
				self._numTypes=p[1]
			if (p[0]=='sequence'):
				self._sequence=p[1]
			if (p[0]=='numChains'):
				self._numChains=p[1]
			if (p[0]=='volFraction'):
				self._volFraction=p[1]
			if (p[0]=='angles'):
				self._angles=bool(p[1])
			if (p[0]=='angle'):
				self._angleString.append(p[1])
				splitString=re.split(':',self._angleString[-1].strip())
				self._angleSequence.append(splitString[0])
				self._angleType.append(int(splitString[1]))
				
	def printDataFilePreview(self):
		f=open(self._nameString,"r")
		for line in f:
			if "Atoms" in line:
				break
			print line,
		f.close()
	
	def printDataFile(self):
		f=open(self._nameString,"r")
		for line in f:
			print line,
		f.close()	
		
				
	def systemName(self):
	#Calculates name used for output data file
		volFractionString=str(self._volFraction).rstrip()
		volFractionString2=volFractionString.replace(".","p")
		self._nameString=str(self._sequence).rstrip()+"_"+str(self._numChains).rstrip()+"_"+volFractionString2+".data"

	def parseSequence(self):
	#Parses sequences. Calculates:
	#	1. nBondTypes - number of different bond types
	#	2. nBondType - List containing bond information for chain (as integer)
	#	3. numBondsPolymer - number of bonds per chain
	#	4. nBondTypesList - original bond sequence from file
	#	5. seqList - List containing sequence
	#	6. seqListN - Number of different types in polymer
	#	7. nTypesinChain - number of monomer types in chain
	
		self._seqList=re.split('-|_|;|:|#',self._sequence)
		self._seqListN=map(int,self._seqList)
		self._nTypesinChain=len(set(self._seqListN))
		self._nBondTypes=' '.join([i for i in self._sequence if not i.isdigit()])
		self._nBondTypesList=re.split(' ',self._nBondTypes.rstrip())
		self._nBondTypes=len(set(self._nBondTypesList))
   
		mapBond=set(self._nBondTypesList)
		mapBond=list(mapBond)
		mapBond2=range(self._nBondTypes+1)
		self._nBondType=[]
   
		#Map bond symbols to bond type numbers
		for x in range(0,len(self._nBondTypesList)):
			ind=mapBond.index(self._nBondTypesList[x])+1
			self._nBondType.append(mapBond2[ind])
	
			self._numBondsPolymer=len(self._nBondType)
			
	def calcSystemParameters(self):
		self._numMonMol=len(self._seqList)
		self._nM=int(self._numMonMol)*int(self._numChains)
		boxVolume=self._nM/(self._density*float(self._volFraction))
		self._boxLength=math.pow(boxVolume,0.33333)
		self._numAtoms=boxVolume*self._density
		self._numWat=self._numAtoms-self._nM
		self._numBonds=float(self._numBondsPolymer)*float(self._numChains)

	def printSystemInfo(self):
		print "Printing Info for system..."
		#print "File %s loaded...\n" %self._fN
		
		print "Input Parameters: "
		print "Number of Atom Types:\t\t%s" %self._numTypes
		print "Polymer sequence:\t\t%s" %self._sequence
		print "Number of polymer chains:\t\t%s" %self._numChains
		print "Polymer volume fraction:\t\t%s" %self._volFraction
		print "System Angles:\t\t%s" %self._angles
		if (self._angles==True):
			print "Number of Angle Types:\t\t%s" %self._nAngleTypes
			#for x in self._angleString:
			#	print x
		
	def initializeDataStructures(self):
		self._atomInformation=np.zeros(shape=(int(self._numAtoms),6))
		self._bondInformation=np.zeros(shape=(int(self._numBonds),4))
		self._angleInformation=np.zeros(shape=(int(self._numAtoms),5))
		
		
	def chainCoordinates(self,ch,startIndex):
	#Format of atom section:
	#atom-ID molecule-ID atom-type x y z
	#startIndex-2 occurs because python data structures start ordering at
	#index 0 and startIndex is iterated after initial assignment
		for x in range(0,(self._numMonMol)):
			self._atomInformation[startIndex-1,0]=startIndex
			startIndex=startIndex+1
			self._atomInformation[startIndex-2,1]=ch+1
			self._atomInformation[startIndex-2,2]=self._seqList[x]
	   
			random.seed()

			if (x==0):
				self._atomInformation[startIndex-2,3]=random.random()*self._boxLength
				self._atomInformation[startIndex-2,4]=random.random()*self._boxLength
				self._atomInformation[startIndex-2,5]=random.random()*self._boxLength
			if (x>0):
				self._atomInformation[startIndex-2,3]=self._atomInformation[startIndex-3,3]+(random.random()-0.5)
				self._atomInformation[startIndex-2,4]=self._atomInformation[startIndex-3,4]+(random.random()-0.5)
				self._atomInformation[startIndex-2,5]=self._atomInformation[startIndex-3,5]+(random.random()-0.5)

	def chainBonds(self,ch,startIndex):
	#Format of bond section:
	#ID type atom1 atom2
		for x in range(0,(self._numBondsPolymer)):
			self._bondInformation[startIndex-1,0]=startIndex
			self._bondInformation[startIndex-1,1]=self._nBondType[x]
			self._bondInformation[startIndex-1,2]=((self._numBondsPolymer+1)*ch)+x+1
			self._bondInformation[startIndex-1,3]=((self._numBondsPolymer+1)*ch)+x+2
			startIndex=startIndex+1
	
	def chainAngles(self,ch,startIndex):
	#Format of angle section:
	#ID type atom1 atom2 atom3
		found=0
		
		for x in range(0,(self._numMonMol-2)):
			interest=self.reconstruct(x)
			for n in range(0,len(self._angleSequence)):
				if interest in self._angleSequence[n]:
					self._angleInformation[startIndex-1,0]=self._numAngles+1
					self._angleInformation[startIndex-1,1]=self._angleType[n]
					self._angleInformation[startIndex-1,2]=((self._numMonMol)*ch)+x+1
					self._angleInformation[startIndex-1,3]=((self._numMonMol)*ch)+x+2
					self._angleInformation[startIndex-1,4]=((self._numMonMol)*ch)+x+3
					startIndex=startIndex+1
					self._numAngles+=1
			
	def reconstruct(self,x):
	  seq=self._seqList[x:x+3]
	  bSeq=self._nBondTypesList[x:x+2]
	  n=""
	  for b in range(0,len(seq)):
		  n+=str(seq[b])
		  if b<2:
			  n+=str(bSeq[b])
	  return n

   
