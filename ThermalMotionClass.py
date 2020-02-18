import math
import numpy as np 


import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.stats import linregress 
from scipy.stats import rayleigh
from scipy.stats import chisquare 

class thermalMotion():
    
    def __init__(self):
        self.xVals = []
        self.yVals = []
        self.pixelToUM = 0.1155
        self.dSqrValsO = []
        self.dt = 0.5
        self.dSqrValsPP = []


    def getRSqr(self,x,y):
        # Computes Distance-Squared on the xy Plane 
        # and rounds to 2 decimal places
        r = (x*x + y*y)
        r = round(r,2)
        return r

    def getFilesData(self,N):
        if N <= 0:
            while True:
                try:
                    N = int(input("How many files would you like to upload?: "))
                    assert(N > 0), 'Number must be bigger than 0:'
                    break
                except:
                    print("Enter a number greater than 0! ")

        fileCount = 1
        while fileCount <= N:
            try:
                fileName = input("Enter the  file name ("+str(fileCount)+" out of %d): "%(N))
                fileObject = open(fileName,'r')
                data = fileObject.readlines()
                for i in range(2,len(data),1):
                    temp = data[i].strip().split("\t")

                    self.xVals += [float(temp[0]) * self.pixelToUM]
                    self.yVals += [float(temp[1]) * self.pixelToUM]
                fileCount += 1
                break
            except OSError:
                print('cannot open', fileName)
        

    def distanceSquaredO(self):
        # -- Calculate Distance-Squared from Initial Position
        initial = [self.xVals[0],self.yVals[0]]
        for j in range(1,len(self.xVals),1):
            dx = self.xVals[j] - initial[0]
            dy = self.yVals[j] - initial[1]
            self.dSqrValsO += [self.getRSqr(dx,dy)]
        return self.dSqrValsO


    def distanceSquaredPP(self):
        # -- Calculates Distance-Squared from Previous Point
        for j in range(1,len(self.xVals),1):
            dx = self.xVals[j] - self.xVals[j-1]
            dy = self.yVals[j] - self.yVals[j-1]
            self.dSqrValsPP += [self.getRSqr(dx,dy)]

        print(self.dSqrValsPP)

    def distSqPlotvsTime(self):
        self.distanceSquaredO()

        timeVals = np.linspace(0,60,len(self.dSqrValsO))
        print(len(timeVals))
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        ax.scatter(timeVals,self.dSqrValsO,color = 'r',marker = 'x',label = "Raw Data")

        # --- Plotting Main Title and Axis Titles --- #

        fig.suptitle("Distance Sqaured vs. Time")
        plt.xlabel('Time (s)')
        plt.ylabel("Distance Squared")
        # --- Getting Fit Data --- # 
        lineInfo = linregress(timeVals,self.dSqrValsO)
        m, b, Rval,stderr = lineInfo[0],lineInfo[1],lineInfo[2],lineInfo[4]
        #print("slope: %.2f, intercept: %.2f, R-val-sq: %.3f, slope err: %.3f" %(m,b,Rval*Rval,stderr))
        # --- Line of Best Fit --- # 
        yLine = timeVals*m + b 
        chiLinear = (chisquare(self.dSqrValsO,yLine))[0]
        redChiLinear = chiLinear/(len(self.xVals) - 1)

        # --- More Plotting --- # 
        plt.plot(timeVals,yLine,'b',label = "Linear Fit: y = " + str(round(m,2)) + "x + " + str(round(b,2)))
        plt.grid(True)
        plt.rc('grid', linestyle=":", linewidth=1, color='gray')
        plt.tick_params(axis='both', which='major', labelsize=10)
        ax.legend(loc = "lower right")
        plt.show()
        
        print("\n")
        # --- Providing User with Fit Data --- # 
        
        print("Fit Data: R-Squared (Correlation Coefficient) = %.3f" %(Rval))
        print("Chi-Squared = %.3f" %(chiLinear))
        print("Reduced Chi-Squared = %.3f" %(redChiLinear))
        print("Slope Error = %.3f" %(stderr))
    
    def rayleigh(self,r,D):
        # -- Rayleigh Distribution -- #
        # D is the diffusion coeffcient
        t = 0.5
        return (r/(2*D*t))*np.exp(-(r*r)/(4*D*t))

    def getHist(self):
        self.distanceSquaredPP()
        # --- data is an array of values such that each element
        # is the distance from its previous point --- #
        
        big = (max(self.dSqrValsPP))
        small = min(self.dSqrValsPP)
        # Making Bin Size
        binSize = big/50
        # --- Initializing bins --- #
        bins = np.linspace(small,big,50, endpoint = True)

        # --- Collecting Frequency -- #
        for i in range(0,len(self.dSqrValsPP),1):
            temp = self.dSqrValsPP[i-1]/binSize
            idx = int(np.floor(temp))
            bins[idx-1] += 1
 
        # --- Plotting the Histogram --- #
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        plt.hist(self.dSqrValsPP,bins = 50, density = True, color = 'red', edgecolor = 'black',label="Histogram of Particle Moving Step Distance - Squared")
        
        plt.title('Particle Step Distances with Time Step: dt = 0.5s ')
        plt.xlabel('Step Distance in Micrometres')
        plt.ylabel('Frequency of Step Distance Obsevered')
        
        

        binSpace = np.linspace(0.0,1.2,50)
        optDVal,dUnc = curve_fit(self.rayleigh,binSpace,bins)
        ax.plot(binSpace,self.rayleigh(binSpace,optDVal), label = "Rayleigh Fit")
        print(dUnc)
        stdDev = np.sqrt(np.diag(dUnc))
        print("standard dev: %.5f" %stdDev)
        plt.legend()
        plt.show()






    
    
