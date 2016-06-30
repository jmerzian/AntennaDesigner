##########################################
#PEEC (Partial Element Equivalent Circuit) based indutor calculator based on:
#http://nvlpubs.nist.gov/nistpubs/jres/69C/jresv69Cn2p127_A1b.pdf
#and 
##########################################

#antenna class, all sizes are in mm
import numpy as nm
import matplotlib.pyplot as plt

#minimum sizes based on manufacturer capabilities
minGap = 0.1
minWidth = 0.1

#the resistivity of the conductive material, assuming copper 1.68×10^-8 ohm/meter
#converted to ohm/mm = 1.68×10^-5
resistivity = 1.68 * nm.power(10,-5,dtype=nm.longdouble)
conductor_thickness = 0.0175

#random utility functions
def num2str(num, precision = 3): 
    return "%0.*f" % (precision, num) 


#trace class, defines what a trace is and how to calculate the characteristics
class Trace:
    def __init__(self, _startPos, _stopPos, _width, _height):
        #physical position of trace
        self.start = _startPos
        self.stop = _stopPos
        #physical paramaters of the trace, converted to cm
        self.width = nm.longdouble(_width/ 10)
        self.length = nm.longdouble(self.GetLength(_startPos,_stopPos)/10)
        self.height = nm.longdouble(_height/ 10)
        #electrical paramaters of the trace
        #calculate the resistance with r*l/(h*w) (ignores thermal effects)
        self.R = nm.longdouble(resistivity * 0.1 * self.length/(self.height*self.width))

        #calculate self inductance
        #.002 l (ln(2l/(w + t))+0.50049+((w+t)/(3l)))
        self.L = nm.longdouble(nm.log(2*self.length/(self.width+self.height)))
        self.L += 0.50049
        self.L += nm.longdouble((self.width+self.height)/(3*self.length))
        self.L *= 0.002*self.length
      
    #get the length of a trace
    def GetLength(self, start, stop):
        return nm.sqrt(nm.square(start[0]-stop[0])+nm.square(start[1]-stop[1]))
       
    #calculate mutual inductance between two traces
    #optional offset paramaters to allow for... offsetting
    def MutualInductance(self,otherTrace, zOffset = 0, xOffset = 0, yOffset = 0):
        #inductance in uH
        L = nm.longdouble(0)      

        #rotate depending on direction vector
        vector = [nm.abs((self.start[0] - self.stop[0])/(10*self.length)),
                   nm.abs((self.start[1] - self.stop[1])/(10*self.length))]
                   
        if(vector[0] == 0):
            x = 1
            y = 0
        else:
            y = 1
            x = 0
        
        #Measurements converted to cm
        #x dimension
        E = ((self.start[x] - (otherTrace.start[x]+xOffset)) + 
            (self.stop[x] - (otherTrace.stop[x]+xOffset)))/20
        a = self.width
        d = otherTrace.width
        
        #y dimension, perpendicular dimension
        l_3 = ((self.start[y] - (otherTrace.start[y]+yOffset)) + 
                (self.stop[y] - (otherTrace.stop[y]+yOffset))) /20
        l_2 = otherTrace.length
        l_1 = self.length
        
        #z dimension
        P = ((self.start[2] - (otherTrace.start[2]+zOffset)) + 
            (self.stop[2] - (otherTrace.stop[2]+zOffset))) /20
        b = self.height
        c = otherTrace.height
        
        #arrays holding the X,Y,Z permutations
        x = [E-a, E+d-a, E+d, E]
        y = [l_3-l_1, l_3+l_2-l_1, l_3+l_2, l_3]
        z = [P-b, P+c-b, P+c, P]
        
        #loop over every permutation of xyz
        for i in range(0,4):
            for j in range(0,4):
                for k in range(0,4):
                    #multiply f(x,y,z) by -1^(1+i+j+k) because... reasons?
                    out = nm.power(-1,i+j+k+2) * self.Mb(x[i],y[j],z[k])
                    L += out
                 
        #scale to the correct values
        L *= nm.divide(0.001,(a*b*c*d))
        return L
        
    def Mb(self,x,y,z):
        #value squares so I don't have to repeatedly compute them
        x2 = nm.square(x)
        y2 = nm.square(y)
        z2 = nm.square(z)
        #value to the third power so I don't have to repeatedly compute them
        x3 = nm.power(x,3)
        y3 = nm.power(y,3)
        z3 = nm.power(z,3)
        #value to the forth power so I don't have to repeatedly compute them
        x4 = nm.power(x,4)
        y4 = nm.power(y,4)
        z4 = nm.power(z,4)
        #distance formula because it's EVERYWHERE
        d = nm.sqrt(x2+y2+z2)
        #breaking rediculously long equation into parts
        #x (y^2 z^2/4 - y^4/24 - z^4/24)
        M_00 = x*(y2*z2/4 - y4/24 - z4/24)
        #M_00 * ln((x+sqrt(x^2+y^2+z^2))/(sqrt(y^2+Z^2)))
        M_01 = nm.log(nm.divide(x+d,nm.sqrt(y2+z2)))
        
        # +
        #y (x^2 z^2/4 - x^4/24 - z^4/24)
        M_10 = y*(x2*z2/4 - x4/24 - z4/24)
        #M_10 * ln((y+sqrt(x^2+y^2+z^2))/(sqrt(x^2+z^2)))
        M_11 = nm.log(nm.divide(y+d,nm.sqrt(x2+z2)))
        
        # +
        #z (x^2 y^2/4 - x^4/24 - y^4/24)
        M_20 = z*(x2*y2/4 - x4/24 - y4/24)
        #M_20 * ln((z+sqrt(x^2+y^2+z^2))/(sqrt(x^2+y^2)))
        M_21 = nm.log(nm.divide(z+d,nm.sqrt(x2+y2)))
            
        # +
        # 1/60 * (x^4 + y^4 + z^4 - 3x^2y^2 - 3y^2z^2-3z^2x^2)*sqrt(x^2+y^2+z^2)
        M_3 = 1/60*(x4 + y4 + z4 - 3*x2*y2 - 3*y2*z2 - 3*z2*x2)*d
        
        #If any two of the variables x, y, and z approach zero, all terms go to zero except the square root term.
        if(nm.isinf(M_01) or nm.isnan(M_01)):
            return M_3
        if(nm.isinf(M_11) or nm.isnan(M_11)):
            return M_3
        if(nm.isinf(M_21) or nm.isnan(M_21)):
            return M_3
            
        # -
        # x*y*z^3/6 * Tan^-1(xy/(z*sqrt(x^2+y^2+z^2)))
        M_4 = x*y*z3/6 * nm.arctan(x*y/(z*d))
        # -
        # x*z*y^3/6 * Tan^-1(xz/(y*sqrt(x^2+y^2+z^2)))
        M_5 = x*z*y3/6 * nm.arctan(x*z/(y*d))
        # -
        # z*y*x^3/6 * Tan^-1(yz/(x*sqrt(x^2+y^2+z^2)))
        M_6 = z*y*x3/6 * nm.arctan(y*z/(x*d))
            
        #If either x or y or z approaches zero, all inverse tangents go to zero.
        
        M = M_00*M_01 + M_10*M_11 + M_20*M_21 + M_3 - M_4 - M_5 - M_6
        
        if(nm.isinf(M) or nm.isnan(M)):
            #print(str(M_00)+"*"+str(M_01))
            #print([x,y,z])
            print("error with mutual inductance calculation...")
            print(str(M_00*M_01)+"+"+str(M_10*M_11)+"+"+str(M_20*M_21)+"+"+str(M_3)+"-"+str(M_4)+"-"+str(M_5)+"-"+str(M_6))
        
        return M
        
#antenna class,
#assumes all traces are connected, IE 0->1->2->3->4 etc.
        #if gap or trace width are not given automaticall space it such that 
        #the width and gap are as wide as possible
        #Currently doesn't work too wel...
class Antenna:       
    def __init__(self, _width, _length, _turns, _gap = -1, _traceWidth = -1):
        #physical paramaters of the antenna
        self.width = _width
        self.length = _length
        self.turns = _turns
        if(_gap == -1):    
            if(_width > _length):
                _gap = (_length/(2*_turns))/2
            else:
                _gap = (_width/(2*_turns))/2
        self.gap = _gap
                
        if(_traceWidth == -1):    
            if(_width > _length):
                _traceWidth = (_length/(2*_turns)) - _gap
            else:
                _traceWidth = (_width/(2*_turns)) -_gap
        
        #physical paramaters of traces
        self.traces = []
        self.trace_Width = _traceWidth
        self.trace_Height = conductor_thickness
        
        #electrical paramaters of the antenna
        self.L = 0
        self.C = 0
        self.R = 0
        self.Z = 0
        
        #special parameters, set outside the main class
        self.Q = 0
        self.layer = 0
        
    def DesignAntenna(self, _layer = 2, thickness = 0.075):
        self.layer = _layer
        for n in range(0,_layer):
            if(self.CoilAntenna(thickness*n,nm.power(-1,n)) == -1):
                return -1
    
    def CoilAntenna(self, z = 0, clockwise = 1):
        #Create all the turns and add them to the self.traces array
        #delta is the total space between two lines
        delta = (self.gap + self.trace_Width)
        if(2*delta*self.turns > self.width or 2*delta*self.turns > self.length):
            print("invalid design")
            return -1
            
        if(clockwise == 1):
            for x in range(0,self.turns):
                #left
                start = [x*delta,x*delta,z]
                stop = [self.width-x*delta,x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #bottom
                start = [self.width-x*delta,x*delta,z]
                stop = [self.width-x*delta,self.length-x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #right
                start = [self.width-x*delta,self.length-x*delta,z]
                stop = [(x+1)*delta,self.length-x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #top
                start = [(x+1)*delta,self.length-x*delta,z]
                stop = [(x+1)*delta,(x+1)*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
        else:
            for x in range(0,self.turns):
                #left
                start = [x*delta,x*delta,z]
                stop = [self.width-x*delta,x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #top
                start = [(x+1)*delta,self.length-x*delta,z]
                stop = [(x+1)*delta,(x+1)*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #right
                start = [self.width-x*delta,self.length-x*delta,z]
                stop = [(x+1)*delta,self.length-x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
                #bottom
                start = [self.width-x*delta,x*delta,z]
                stop = [self.width-x*delta,self.length-x*delta,z]
                newTrace = Trace(start, stop, self.trace_Width, self.trace_Height)
                #add the trace so the antenna characteristics can be calculated
                self.traces.append(newTrace)
                
        #add each trace parameter to coil parameter
        for trace in self.traces:
            self.L += trace.L
            self.R += trace.R
        
        #which trace the iterator is currently on
        x = 0
        traceCount = len(self.traces)
        Mplus = 0
        Mminus = 0
        for i in range(0,traceCount):
            #increment x to prevent double counting traces    
            x += 1
            for j in range(x,traceCount):
                #make sure the traces are paralell
                if(i%2 == j%2):
                    #check they are going the same direction
                    L = self.traces[i].MutualInductance(self.traces[j])
                    if(i%4 == j%4):
                        Mplus += L
                    else:
                        Mminus += L      
        
        self.L += Mplus - Mminus

    #TODO this might not actually be correct... however I don't have a good way to test ATM
    #Seems pretty good, used in other calculations and everything looks okay...
    #optional offset values to allow for... offsetting
    def Mutual(self, otherAntenna, zOffset = 1, xOffset = 0, yOffset = 0):
        Mplus = 0
        Mminus = 0        
        
        traceCount = len(self.traces)
        otherTraces = len(otherAntenna.traces)
        for i in range(0,traceCount):
            #couple with all the traces of the other antenna
            for j in range(0,otherTraces):
                #make sure the traces are paralell
                if(i%2 == j%2):
                    #Get the mutual inductance between two traces
                    L = self.traces[i].MutualInductance(otherAntenna.traces[j],
                                                        zOffset,xOffset,yOffset)
                    #Check they are going the same direction, makes some bad assumptions...                    
                    if(i%4 != j%4):
                        Mplus += L
                    else:
                        Mminus += L 
                      
        #print(str(Mplus) + " - " + str(Mminus))
        return Mplus + Mminus
    
    #k = M/(sqrt(L1*L2))
    def K(self, otherAntenna, zOffset = 1, xOffset = 0, yOffset = 0):
        Mutual = self.Mutual(otherAntenna,zOffset,xOffset,yOffset)
        return Mutual/(nm.sqrt(self.L*otherAntenna.L))
        
    def GetDimensions(self):
        return [self.width,self.length,self.turns,self.gap,self.trace_Width]
    
    def Draw(self, figure, label=-1):
        plt.figure(figure) 
        if(label != -1):
            plt.title(label)
        else:
            plt.title("trace:"+str(self.trace_Width) + "x" + str(self.gap)+"antenna:" + str(self.width) + "x" + str(self.length) + "x" + str(self.turns))
        plt.axes()
        for i in range(0,len(self.traces)):
            vector = [(self.traces[i].start[0] - self.traces[i].stop[0])/(self.traces[i].length*10),
                      (self.traces[i].start[1] - self.traces[i].stop[1])/(self.traces[i].length*10)]
                   
            if(vector[0] == 0):
                rectangle = plt.Rectangle(self.traces[i].start,
                                          self.traces[i].width*10,
                                          self.traces[i].length*10*-vector[1], fc='r')
            else:
                rectangle = plt.Rectangle(self.traces[i].start,
                                          self.traces[i].length*10*-vector[0],
                                          self.traces[i].width*10, fc='r')
            plt.gca().add_patch(rectangle)
            #old line method...
            #line = plt.Line2D([x1,x2],[y1,y2],lw=2.5)
            #plt.gca().add_line(line)
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.axis("scaled")
        plt.show()
        
    def GenerateEagle(self,name = "Test"):    
        #Eagle wire command WIRE CW (0 1) (0 -1);        
        output = open("Output/"+str(name)+".scr","w")
        output.write("GRID MM \n")
        output.write("LAYER 1 \n")
        
        #print out all the reaces
        for trace in self.traces:
            output.write("WIRE " + num2str(trace.width*10) + " 'antenna' ( "+
                           num2str(trace.start[0])+" "+num2str(trace.start[1])+
                        ") ("+num2str(trace.stop[0])+" "+num2str(trace.stop[1])+" ); \n")

        output.close()
        
    def GenerateRoundEagle(self,name = "Test"):    
        #Eagle wire command WIRE CW (0 1) (0 -1);   
        #Eagle arc command ARC CW (0 1) (0 -1) (1 0);
        output = open("Output/"+str(name)+".scr","w")
        output.write("GRID MM \n")
        output.write("LAYER 1 \n")
        
        #iterates over the traces,
        #TODO make this better, assuming a standard 4 turn antenna...
        delta = self.trace_Width + self.gap
        
        #print out all the reaces
        for n in range(0,self.turns):
            #height of the arc
            arc = nm.abs(((delta*n) - (self.width - delta*n))/2)
            
            #left
            x0 = num2str(delta*n)
            x1 = num2str(delta*n)
            y0 = num2str(delta*(n-1))
            y1 = num2str(self.length - delta*n - arc)
            output.write("WIRE " + num2str(self.trace_Width) + 
                " '"+name+"' ( " + x0 + "  " + y0 + " )" +
                "( " + x1 + "  " + y1 + " ); \n")
                
            
            #top left arc
            x0 = num2str(self.width - delta*n)
            x1 = num2str(delta*n)
            y0 = num2str(self.length - delta*n  - arc)
            y1 = num2str(self.length - delta*n  - arc)
            output.write("ARC '"+name+"' CCW (" +
                    x0 + " " + y0 + ") (" +
                    x1 + " " + y1 + ") (" +
                    x1 + " " + y1 + "); \n")
                    #num2str(self.length - delta*n + arc) + "); \n")
            
               
            #right
            x0 = num2str(self.width - delta*n)
            x1 = num2str(self.width - delta*n)
            y0 = num2str(delta*n)
            y1 = num2str(self.length - delta*n - arc)
            output.write("WIRE " + num2str(self.trace_Width) + 
                " '"+name+"' ( " + x0 + " " + y0 + " )"+
                "( " + x1 + " " + y1 + " ); \n")
                
            #bottom
            x0 = num2str(self.width - delta*n)
            x1 = num2str(delta*(n+1))
            y0 = num2str(delta*n)
            y1 = num2str(delta*n)
            output.write("WIRE " + num2str(self.trace_Width) + 
                " '"+name+"' ( " + x0 + " " + y0 + " )"+
                "( " + x1 + " " + y1 + " ); \n")
        
        #if it's dual layer draw the second layer, mirrored
        if(self.layer > 1):
            #define bottom layer
            output.write("LAYER 16 \n")
            #print out all the traces
            for n in range(0,self.turns):
                #height of the arc
                arc = nm.abs(((delta*n) - (self.width - delta*n))/2)
                
                #left
                x0 = num2str(delta*n)
                x1 = num2str(delta*n)
                y0 = num2str(delta*n)
                y1 = num2str(self.length - delta*n - arc)
                output.write("WIRE " + num2str(self.trace_Width) + 
                    " '"+name+"' ( " + x0 + "  " + y0 + " )" +
                    "( " + x1 + "  " + y1 + " ); \n")
                    
                
                #top left arc
                x0 = num2str(self.width - delta*n)
                x1 = num2str(delta*n)
                y0 = num2str(self.length - delta*n  - arc)
                y1 = num2str(self.length - delta*n  - arc)
                output.write("ARC '"+name+"' CCW (" +
                        x0 + " " + y0 + ") (" +
                        x1 + " " + y1 + ") (" +
                        x1 + " " + y1 + "); \n")
                        #num2str(self.length - delta*n + arc) + "); \n")
                
                   
                #right
                x0 = num2str(self.width - delta*n)
                x1 = num2str(self.width - delta*n)
                y0 = num2str(delta*(n-1))
                y1 = num2str(self.length - delta*n - arc)
                output.write("WIRE " + num2str(self.trace_Width) + 
                    " '"+name+"' ( " + x0 + " " + y0 + " )"+
                    "( " + x1 + " " + y1 + " ); \n")
                    
                #bottom
                x0 = num2str(self.width - delta*(n+1))
                x1 = num2str(delta*n)
                y0 = num2str(delta*n)
                y1 = num2str(delta*n)
                output.write("WIRE " + num2str(self.trace_Width) + 
                    " '"+name+"' ( " + x0 + " " + y0 + " )"+
                    "( " + x1 + " " + y1 + " ); \n")
        
        output.close()
        
     
"""
def Test():
    
    mainAnt = Antenna(30,40,2,0.3,0.5)  
    mainAnt.DesignAntenna()
    k = []
    #code to test the antenna coupling calculation
    for n in range(1,5):
        testAnt = Antenna(30,40,n,0.3,0.5)
        testAnt.DesignAntenna()
        #testAnt.DrawAnt()
        #print("N: " + str(n) + " L: " + str(testAnt.L) + " R: " + str(testAnt.R))
        ktemp = []
        for d in range(1,10):
            M = mainAnt.Mutual(testAnt,d)
            K = M/(nm.sqrt(testAnt.L*mainAnt.L))
            #print(str(K) + " = " + str(M) + "/ sqrt" + str(testAnt.L) + " " + str(mainAnt.L))
            ktemp.append(K)
        k.append(ktemp)
    
    plt.figure(1)    
    
    plt.subplot(221)
    plt.title("test case 1")
    for n in range(0,len(k)):
        plt.plot(k[n])
    plt.grid(True)
    
    plt.show()
    
Test()
"""
"""
def Test():
    #Code to test the offset of the antenna calculations
    #antenna values
    ant0 = []
    #30mm x 40mm, 0-10 turns, 0.3mm gap, 0.5mm width
    #Panasonic
    testData0 = [0.143,0.451,0.872,1.38,1.956,2.583,3.250,3.943,4.654,5.372]
    #ST
    testData1 = [0.13328,0.42085,0.81101,1.28,1.81,2.36,2.91,3.47,4.02,4.57]
    for n in range(1,11):
        testAnt = Antenna(30,40,n,0.3,0.5)
        if(testAnt.DualLayer() == -1):
            break
        #testAnt.DrawAnt()
        print("N: " + str(n) + " L: " + str(testAnt.L) + " R: " + str(testAnt.R))
        #add the actual inductance
        ant0.append(testAnt.L)
    
    #antenna values
    ant1 = []
    #30mm x 40mm, 0-10 turns, 0.3mm gap, 0.5mm width
    #Panasonic
    testData2 = [0.074,0.230,0.435,0.671,0.926,1.189,1.451,1.702,1.937,2.107]
    #ST
    testData3 = [0.06666,0.19876,0.36427,0.54412,0.72448,0.89215,1.04,1.16,1.26,1.33]
    for n in range(1,11):
        testAnt = Antenna(20,20,n,0.3,0.5)
        if(testAnt.DualLayer() == -1):
            break
        testAnt.Draw("antenna"+str(n))
        print("N: " + str(n) + " L: " + str(testAnt.L) + " R: " + str(testAnt.R))
        #add the actual inductance so it can be compared
        ant1.append(testAnt.L)
        
    plt.figure(1)
    #the raw data inputs
    err0 = nm.array([testData0,testData1,testData1])
    err1 = nm.array([testData2,testData3,testData3])
    
    #the mean value
    mean0 = nm.mean(err0, axis = 0)
    mean1 = nm.mean(err1, axis = 0)
    
    #the mean averaged whatever thingy
    err0 = nm.divide(err0,mean0)
    err1 = nm.divide(err1,mean1)
    
    plt.subplot(221)
    plt.title("test case 1")
    plt.plot(ant0)
    plt.plot(testData0)
    plt.plot(testData1)
    plt.grid(True)
    
    plt.subplot(222)
    plt.title("test case 2")
    plt.plot(ant1)
    plt.plot(testData2)
    plt.plot(testData3)
    plt.grid(True)
    
    plt.subplot(223)
    plt.title("test case 1 mean error")
    plt.plot(stats.sem(err0))
    plt.grid(True)
    
    plt.subplot(224)
    plt.title("test case 2 mean error")
    plt.plot(stats.sem(err1))
    plt.grid(True)
    
    plt.show()    

Test()
"""