from Antenna import Antenna#custom antenna function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as nm
import time
import re

##########################################
#NFC design calculator based on
#http://www.jocm.us/uploadfile/2013/0927/20130927043742253.pdf
##########################################
#Natural frequency in MHz
w_0 = 2 * nm.pi * 13.56
#Equal to V^2/P
R_L = 1250

#next open figure
nFig = 5

#data container for graph stuffs
Q_2 = []
N_2 = []
M_12 = []
V_2 = []

#Get the reccomended turn count for a given geometry
def GetN(L,R, R_l = R_L):
    turns = nm.power(2*R*R_l/(nm.square(w_0)*nm.square(L)),1/3)
    return nm.ceil(turns)

#Get the quality facotr of a coupled tag
#Q_2 = 1/(w_0*L_2/R_l + R_2/(w_0*L_2))
def GetQ(L,R,R_l = R_L):
    Q = 1/((w_0*L/R_l) + (R/(w_0*L)))
    return Q

#Get the power approximation
#R_t = w_0*k^2*L_1*Q_2
def GetR_t(k,L,Q):
    R_t = w_0*nm.square(k)*L*Q
    return R_t
    
#Get the root-mean-square power dissipated by the tag
#P = I^2*R
def GetP(R_t,I = 0.06):
    P = nm.square(I)*R_t
    return P
    
#W_0*M*I*Q
def GetV(R_t,I = 0.06):
    V = I*R_t
    return V

#Find the best antenna design for a given reader and available dimension
#returns the best antenna design
def GetBest(readAnt,width,length,zOffset = 20, plot = 0):
    #hold the R_t map
    R_tMap = []
    
    #Position of max energy for antenna
    max_Pos = []
    #highest energy for an  antenna
    R_tmax = []
    #average energy within the space
    avgR_t = []
    #Consistency of field
    R_tDev = []
    
    bestAnt = []
    #loop through different turn settings
    for n in range(1,4):
        #create a new antenna with the given parameters
        testAnt = Antenna(width, length, n, 0.15)
        testAnt.DesignAntenna()
        #print out the calculated values of the created antenna
        print("N: " + str(n) + " L: " + str(testAnt.L) + " R: " + str(testAnt.R))
        
        #Get the quality factor
        Q = GetQ(testAnt.L,testAnt.R)
        #add it to the antenna if needed later
        testAnt.Q = Q
        
        #Simulation box
        #Currently set so the tag is placed on every corner of the reader
        minXY = [-readAnt.width-width,-readAnt.length-length] 
        #minXY = [-width,-length]
        maxXY = [2*(readAnt.width+width),2*(readAnt.length+length)] 
        #maxXY = [readAnt.width+width,readAnt.length+length]
        
        #Map the R_t of the antenna accross the read range, check every cm
        #readAnt,testAnt,minXY,maxXY,step = 10,figure = 1,zOffset = 20
        #figure empty for no graph
        R_tTemp = offsetMap(readAnt,testAnt,minXY,maxXY,5,5*int(length/10)+n)
        R_tMap.append(R_tTemp)
        
        #Position of min/max energy for antenna
        max_Pos.append(nm.argmax(R_tTemp))
        #lowest/highest energy for an  antenna
        R_tmax.append(nm.amax(R_tTemp))
        #average energy within the space
        avgR_t.append(nm.mean(R_tTemp))
        #Consistency of field
        R_tDev.append(nm.std(R_tTemp))
        
        #Add the antenna to the running for best antenna
        bestAnt.append(testAnt)
    
    if(plot == 1):        
        plt.figure(1)
        
        #Antenna impedance
        plt.subplot(221)
        plt.title("R_t Max")
        plt.plot(R_tmax,label = str(width)+"x"+str(length))
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #Antenna impedance
        plt.subplot(222)
        plt.title("R_t Average")
        plt.plot(avgR_t,label = str(width)+"x"+str(length))
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
        #Antenna impedance
        plt.subplot(223)
        plt.title("Standard Deviation")
        plt.plot(R_tDev,label = str(width)+"x"+str(length))
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
    #WHO WON, WHOS NEXT... I have no idea how to decide...
    print("antenna " + str(nm.argmax(R_tmax)) + "WINS!!!")
    return bestAnt[nm.argmax(R_tmax)]
    
def GetQKRN(readAnt,testAnt,xOffset = 0,yOffset = 0,zOffset = 20, layers = 1):
    #Create the antenna, if it's invalid stop!   
    if(testAnt.DesignAntenna(layers) == -1):
        return [0,0,0,0]
    
    #Get the quality factor and add it to the antenna if needed later
    Q = GetQ(testAnt.L,testAnt.R)
    testAnt.Q = Q
    
    #Get the coupling coeffecient at zOffset centered on the reader
    K = nm.abs(readAnt.K(testAnt,zOffset,xOffset,yOffset))
    #Get the R_t
    R_t = GetR_t(K,readAnt.L,testAnt.Q)
    #get the impedance and add it
    N = GetN(testAnt.L,testAnt.R)
    
    print([Q,K,R_t,N])
    
    return [Q,K,R_t,N]
        
#QKR label for the dataset, figure to add data to
def QKRNGraph(Q_2,k,R_t,N,title = "Q K R",fig = 2):
    plt.figure(fig) 
    
    #Q factor of antenna
    plt.subplot(221)
    plt.title("Q")
    plt.plot(Q_2,label = title)
    plt.xlabel(fig)
    plt.grid(True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    #coupling between A and B
    plt.subplot(222)
    plt.title("k")
    plt.plot(k,label = title)
    plt.xlabel(fig)
    plt.grid(True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    #Antenna impedance
    plt.subplot(223)
    plt.title("R_t")
    plt.plot(R_t,label = title)
    plt.xlabel(fig)
    plt.grid(True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
    
    #ideal N
    plt.subplot(224)
    plt.title("N")
    plt.plot(N,label = title)
    plt.xlabel(fig)
    plt.grid(True)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
    
    
#Creates a 3d array mapping the impedance at points above the antenna
#Returns the heightmap 
def offsetMap(readAnt,testAnt,minXY,maxXY,step = 10,figure = -1,zOffset = 20):
    #X,Y,Z coordinates for the map
    x1 = nm.arange(minXY[0],maxXY[0],step)
    y1 = nm.arange(minXY[1],maxXY[1],step)
    
    X, Y = nm.meshgrid(x1,y1)
    #Calculate the R_t at a given (X,Y) location
    def F(x,y):
        #Get the coupling coeffecient at zOffset centered on the reader
        tempk = nm.abs(readAnt.K(testAnt,zOffset,x,y))
        #Get the R_t
        R_t = GetR_t(tempk,readAnt.L,testAnt.Q)
        print(str([x,y])+"="+str(R_t))
        return R_t  
        
    #loop through possible positions with a step size of step
    z = nm.array([F(x,y) for x,y in zip(nm.ravel(X), nm.ravel(Y))])
    Z = z.reshape(X.shape)
     
    #save the data to a files
    nm.savez("Output/R_tMap"+str(figure),X,Y,Z)
    if(figure != -1):
        Plot3D(X,Y,Z,figure)
    return Z
    
#Plots an XYZ file based on an input, assumes it's a numpy files
def Plot3DFile(file, fig = 1):
    #Get the input file
    data = nm.load("Output/"+str(file)+".npz")
    X = data['X']
    Y = data['Y']
    Z = data['Z']
    data.close()
    #Plot the resultant data
    Plot3D(X,Y,Z,fig)
    
#plot a 3d X Y Z graph on figure fig
def Plot3D(X,Y,Z,fig=1):
    #create a new 3d graph
    fig = plt.figure(fig)
    ax = fig.gca(projection='3d')
    
    #set the graph to be a 3d surface
    surf = ax.plot_surface(X,Y,Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    
    #add a color map to the contour as Z increases
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    #put a color bar on the side to help with determining values
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
def SavePlotToFile(filename,figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.set_size_inches(18.5, 10.5, forward=True)
        fig.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    
#iterate over the given parameters to find the best designs
#Finds best designs for every intermediarry as well
def Iterate(width = [10,45,5],length = [8,9,1],
            turns = [1,11],traceWidth = [0.15,1.65,0.1],
            layers = 1):
     #create a model of the reader antenna
    readAnt = Antenna(80,60,4,0.3,1)
    readAnt.DesignAntenna()
    #readAnt.Draw(0)
    
    #Holds the best antenna design
    TheBest = []
    #Create an output file for the best designs
    bestFile = open('TheBest.txt', 'w')
    
    antennaSave = []
    
    #number of layers (1 or 2)
    for L in range(1,layers+1):
        #length of the antennas
        for l in nm.arange(length[0],length[1],length[2]):
            #width of the antenna
            for w in nm.arange(width[0],width[1],width[2]):
                bestR = -1
                bestK = -1
                #temporary placeholder antennas
                bestR_ant = Antenna(80,60,4,0.3,1)
                bestK_ant = Antenna(80,60,4,0.3,1)
                #thickness of traces
                for t in nm.arange(traceWidth[0],traceWidth[1],traceWidth[2]): 
                    #temporary holder for QKR
                    Q = []
                    K = []
                    R = []
                    N = []
                    
                    print(str(l)+"x"+str(w)+"w x" + str(L) +" "+ "g:0.15"+"th"+str(t))
                    #number of turns
                    for n in range(turns[0],turns[1]):
                        ant1 = Antenna(l,w,n,0.15,t)
                        #dual layer, the number of turns is 2*n...
                        QKRN = GetQKRN(readAnt,ant1,layers = L,
                                       xOffset = 40,yOffset = 30)
                        Q.append(QKRN[0])
                        K.append(QKRN[1])
                        R.append(QKRN[2])
                        N.append(QKRN[3])
                        #check if it has the highest K
                        if(QKRN[1]>bestK):
                            #if so make it the new contender and see if anything else can beat it
                            bestK = QKRN[1]
                            bestK_ant = ant1
                            
                        #check if it has the highest R
                        if(QKRN[2]>bestR):
                            #if so make it the new contender and see if anything else can beat it
                            bestR = QKRN[2]
                            bestR_ant = ant1
            
                        #QKR label for the dataset, figure to add data to
                        QKRNGraph(Q,K,R,N,"th:"+str(t),str(l)+"x"+str(w)+"th-"+str(t))
                        #Add the data to the master save file so you don't have to 
                        #simulate every GOD DAMN TIME!!!!
                        antennaSave.append([ant1.GetDimensions(),Q,K,R,N])
    
                SavePlotToFile("Output/"+str(l)+"x"+str(w)+" x " + str(L) + "_40x_30y.pdf")
                
                #add the new best to... the best
                TheBest.append(bestR_ant)
                TheBest.append(bestK_ant)
            
    #Save all the numbers from the trial run
    nm.savez("Output/DataDump"+time.strftime("%d_%m_%Y-%H_%M_%S"),antennaSave)
        
    i = 0
    #save the best designs
    for Best in TheBest:
        bestFile.write("w="+str(Best.width)+" l="+str(Best.length)+" n="+
            str(Best.turns)+" g="+str(Best.gap)+" w="+str(Best.trace_Width)+
            " L="+str(Best.layer)+" R="+str(Best.R)+"\n")
        Best.GenerateRoundEagle("scr/Best"+str(Best.width)+"x"+str(Best.length)+"x"+str(Best.layer)+"_"+str(i))
        i += 1
    bestFile.close()
        
    """
    #Currently set so the tag is placed on every corner of the reader
    minXY = [-readAnt.width,-readAnt.length] 
    #minXY = [-width,-length]
    maxXY = [2*readAnt.width,2*readAnt.length] 
    
    f = 0
    
    for Best in TheBest:
        for z in range(1,101,20):
            name = ("w="+str(Best.width)+" l="+str(Best.length)+" n="+
                str(Best.turns)+" g="+str(Best.gap)+" w="+str(Best.trace_Width)+
                " L="+"{0:.2f}".format(Best.L)+" R="+"{0:.2f}".format(Best.R)+"\n")
            #(readAnt,testAnt,minXY,maxXY,step = 10,figure = -1,zOffset = 20)
            offsetMap(readAnt,Best,minXY,maxXY,10,name,z)
        #save the plots to a single file
        SavePlotToFile("Output/3d" + str(f) + ".pdf")
        f += 1
    """
    
#Main function
def __Main__():
    #iterate over a width,length,turns,trace width,layers
    #Iterate([50,100,10],[50,100,10],[1,15],[0.15,2,0.05],2)
    """
    bestFile = open('AntennaDesigns.txt', 'r')
    iterator = 0
    for line in bestFile:
        string = re.sub("[a-zA-Z]+(=)", " ", line)
        x = [float(i) for i in string.split()]
        Best = Antenna(x[0],x[1],int(x[2]),x[3],x[4])
        Best.DesignAntenna(int(x[5]))
        name = "scr/Best"+str(Best.width)+"x"+str(Best.length)+"x"+str(Best.layer)+"_"+str(iterator)
        Best.GenerateRoundEagle(name)
        iterator += 1
        print(x)
    
    """
    #_width, _length, _turns, _gap = -1, _traceWidth = -1
    antenna = Antenna(6,40,2,0.15,1.15)
    antenna.DesignAntenna(2)
    antenna.GenerateRoundEagle("RoundedAntenna")
       
    plt.show()
    
__Main__()