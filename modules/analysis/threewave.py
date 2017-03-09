#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
# CoCoPy - A python toolkit for rotational spectroscopy
#
# Copyright (c) 2013 by David Schmitz (david.schmitz@chasquiwan.de).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the “Software”), to deal in the 
# Software without restriction, including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# MIT Licence (http://mit-license.org/)
#
################################################################################
#
# Module, threewave.py, written by V. Alvin Shubert, (alvinshubert@gmail.com)
#
################################################################################

import numpy as np
import spfit

'''
Todo:
    - 1 
    - 2 
    - 3 
'''
class threewave:
    '''
    The threewave class tries to cover everything, which is related to searching for
    three-wave-mixing cycles.
    '''
    
    def __init__(self, fname = 'unnamed.txt', fegy = 'unnamed'):
        '''
        Starts an empty threewave instance.

        Parameters
        ----------
        fname (str): Filename of three-wave mixing file ('unnamed.txt').
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        self.hello = 'hello'
        self.fname = fname
        self.fegy = fegy
        self.cycles = np.array([])
        self.egy_mhz = np.array([])
        self.egy_qn = np.array([])
        self.twists = np.array([])
        self.listens = np.array([])
        self.drives = np.array([])
                
    def save(self, fname=''):
        '''
        Saves the cycles in a textfile.
        The quantum numbers, transition frequencies, and log(intensities) are
        saved into the a textfile named fname in the format:
        
        J1, ka1, kc1, J2, ka2, kc2, J3, ka3, kc3, frq1, frq2, frq3, int1, int2, int3
        
        where: 
        
        (J1, ka1, kc1) -> (J2, ka2, kc2), frq1, int1 are for drive parameters
        (J2, ka2, kc2) -> (J3, ka3, kc3), frq2, int2 are the twist parameters
        (J3, ka3, kc3) -> (J1, ka1, kc1), frq3, int3 are the listen parameters

        Parameters
        ----------
        fname (str): Filename of cycles file ('').
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if fname != '':
            self.fname = fname
            
        if (len(self.drives)<2):
            print 'Sorry, nothing to save. Try searching for some cycles, maybe.'
            return
        #output = np.column_stack((self.cycles))
        
        dtl=np.arange(15)*1.0
        for x in np.arange(len(self.drives)):
            d = self.drives[x]
            t = self.twists[x]
            l = self.listens[x]
            
            dtl[0:6] = d[0:6]
            dtl[9] = d[6]
            dtl[10] = t[6]
            dtl[11] = l[6]
            dtl[12] = d[7]
            dtl[13] = t[7]
            dtl[14] = l[7]
            if all(d[0:3]==t[0:3]):
                dtl[6:9]=t[3:6]
            else:
                dtl[6:9]=t[0:3]
            if x == 0:
                output = dtl
            else:
                output = np.vstack((output,dtl))
        
        np.savetxt(self.fname, output, fmt='%3d,%3d,%3d,%3d,%3d,%3d,%3d,%3d,%3d'
                                            ',%10.4f,%10.4f,%10.4f,%8.4f,%8.4f,%8.4f')


    def savepdf(self, fname=''):
        '''
        Saves the cycles into a tex file that is used to generate a pdf. 
        Transition type are color coded and the predicted most intense three-
        wave mixing cycle is italicized.
        
        A .txt file is also produced that contains the cycles.
        
        Parameters
        ----------
        fname (str): Filename of cycles file ('').
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        import subprocess

        if fname != '':
            self.fname = fname
            
        if (len(self.drives)<2):
            print 'Sorry, nothing to save. Try searching for some cycles, maybe.'
            return

        # Search for the most intense cycle (i.e. such that sum(log(intensities))
        # is a maximum.
        mint= -100.
        for x in np.arange(len(self.twists)):
            if (self.twists[x][7]+self.drives[x][7]+self.drives[x][7] > mint):
                mint = self.twists[x][7]+self.drives[x][7]+self.drives[x][7]
                mintx = x

        # Create the name of the .tex file.
        ftex = self.fname[0:len(self.fname)-3]+'tex'
        
        # The name of the spcat simulation is written to the pdf. Since it might
        # contain special characters, like '_', this section is here to add the
        # necessary '\' for use in LaTex.
        spcat_fname=''
        for x in self.fegy:
            if x != '_':
                spcat_fname = spcat_fname+x
            else:
                spcat_fname = spcat_fname+'\\_'
        
        # Open the LaTex file for writing and write the preamble.
        fp = open(ftex,'w')
        fp.write('\\documentclass[a4paper,english,latex,twoside]{article}\n')
        fp.write('\\newcommand{\\tab}[1]{\\hspace{.05\\textwidth}\\rlap{#1}}\n')
        fp.write('\\usepackage{babel}\n')
        fp.write('\\usepackage[usenames,dvipsnames]{color}')

        # Begin the main part of the LaTex document.Here a title is written, along
        # with an explanation of the file contents and the search limits.
        fp.write('\\begin{document}\n')
        fp.write('\\begin{center}\\begin{Large}\n \\textbf{Three-wave-mixing ' 
                'Cycles for $'+spcat_fname+'$}\\end{Large}\\end{center}\n')
        fp.write('This file presents the three-wave-mixing cycles from the '
                'spcat prediction of $'+spcat_fname+'$. The log(intensity) for each '
                'transition is given in ().')
        fp.write('The transitions are color coded as {\\color{blue}a-type}, '
                '{\\color{OliveGreen}b-type}, and {\\color{Fuchsia}c-type}. '
                'The cycle with the highest combined intensity')        
        fp.write(' is shown in italics.\n\n')        
        fp.write('The search limits were:\n\n')
        fp.write('\\vspace*{.5\\baselineskip}')
        fp.write('   '+str(self.Jmin)+' $\leq J \leq$ '+str(self.Jmax)+'\n\n')
        fp.write('   '+str(self.rfmin)+' MHz $\leq\\nu_{twist}\leq$ '+str(self.rfmax)+' MHz\n\n')
        fp.write('   '+str(self.dmin)+' MHz $\leq\\nu_{drive}\leq$ '+str(self.dmax)+' MHz\n\n')        
        fp.write('   '+str(self.lmin)+' MHz $\leq\\nu_{listen}\leq$ '+str(self.lmax)+' MHz\n \n \n')    
        fp.write('\\vspace*{.5\\baselineskip}')

        # Open the text file for writing. 
        f = open(self.fname,'w')
        f.write('Cycles from the spcat prediction of '+self.fegy+'.\n')
        f.write('The search limits were:\n')
        f.write('   '+str(self.Jmin)+' >= J <= '+str(self.Jmax)+'\n')
        f.write('   '+str(self.rfmin)+' MHz >= twist transition <= '+str(self.rfmax)+' MHz\n')
        f.write('   '+str(self.dmin)+' MHz >= drive transition <= '+str(self.dmax)+' MHz\n')        
        f.write('   '+str(self.lmin)+' MHz >= listen transition <= '+str(self.lmax)+' MHz\n \n')    
        f.write('Below, the log(intensity) for each transition is given in parentheses.\n \n')


        # The next section writes the cycles to the file. The cycles are organized
        # so that all with the same twist transition are written together. Thus
        # the array tp contains the parameters of the previously written twist.
        # If the current twist is the same as the previous, then just the drive and 
        # listens are written. Otherwise a new twist section is started by writing
        # the twist, followed by the drives and listens.
        tp = np.arange(8)*1.0
        for x in np.arange(len(self.twists)):
            d = self.drives[x]
            t = self.twists[x]
            l = self.listens[x]
            # This section determines the transition types so that they can 
            # be appropriately color-coded.
            dt = self.checktranstype(d[0],d[1],d[2],d[3],d[4],d[5])
            tt = self.checktranstype(t[0],t[1],t[2],t[3],t[4],t[5])
            lt = self.checktranstype(l[0],l[1],l[2],l[3],l[4],l[5])
            ttt = ['blue','OliveGreen','Fuchsia']
            tr = ['a','b','c']
            tca = ['black','black','black']
            tra = [dt,tt,lt]
            for y in np.arange(len(tra)):
                for z in np.arange(len(tr)):
                    if tra[y] == tr[z]:
                        tca[y] = ttt[z]
            dc = tca[0]
            tc = tca[1]
            lc = tca[2]
            
            # This section creates strings from the quantum numbers in the
            # found cycles. This procedure is done for easier handling of the
            # information when writing to the files.
            tq = [str(int(y)) for y in t]
            dq = [str(int(y)) for y in d]
            lq = [str(int(y)) for y in l]
            tlow = tq[0]+' '+tq[1]+' '+tq[2]
            tlop = tq[0]+'_{'+tq[1]+tq[2]+'}'
            thgh = tq[3]+' '+tq[4]+' '+tq[5]
            thgp = tq[3]+'_{'+tq[4]+tq[5]+'}'
            dlow = dq[0]+' '+dq[1]+' '+dq[2]
            dlop = dq[0]+'_{'+dq[1]+dq[2]+'}'
            dhgh = dq[3]+' '+dq[4]+' '+dq[5]
            dhgp = dq[3]+'_{'+dq[4]+dq[5]+'}'
            llow = lq[0]+' '+lq[1]+' '+lq[2]
            llop = lq[0]+'_{'+lq[1]+lq[2]+'}'
            lhgh = lq[3]+' '+lq[4]+' '+lq[5]
            lhgp = lq[3]+'_{'+lq[4]+lq[5]+'}'
            
            # Write the cycles to the files. Note f is the text file object and
            # fp is the LaTex file object.
            if any(tp != t):
                f.write('Twist: '+tlow+' <-> '+thgh+' , '+('{0:8.3f}'.format(abs(t[6])))
                        +' MHz, ('+('{0:6.4f}'.format(t[7]))+')\n\n' )
                fp.write('\\vspace*{.5\\baselineskip}')
                fp.write(' \\noindent {\\color{'+tc+'}Twist: $'+tlop+' \leftrightarrow '
                        +thgp+'$ , '+('{0:8.3f}'.format(abs(t[6])))+' MHz, ('
                        +('{0:6.4f}'.format(t[7]))+')}\n\n' )
            f.write('     Drive: '+dlow+' <-> '+dhgh+' , '+('{0:8.3f}'.format(abs(d[6])))
                    +' MHz, ('+('{0:6.4f}'.format(d[7]))+')')
            if x == mintx:
                fp.write('\\textit{')
            fp.write('\\noindent \\tab{{\\color{'+dc+'}Drive: $'+dlop+'  \leftrightarrow  '
                    +dhgp+'$ , '+('{0:8.3f}'.format(abs(d[6])))+' MHz, ('
                    +('{0:6.4f}'.format(d[7]))+')}')
            f.write(' Listen: '+llow+' <-> '+lhgh+' , '+('{0:8.3f}'.format(abs(l[6])))
                    +' MHz, ('+('{0:6.4f}'.format(l[7]))+')\n\n')
            fp.write(' {\\color{'+lc+'}Listen: $'+llop+' \leftrightarrow '+lhgp+'$ , '
                    +('{0:8.3f}'.format(abs(l[6])))+' MHz, ('+('{0:6.4f}'.format(l[7]))+')}')
            if x == mintx:
                fp.write('}')
            fp.write('}\n\n')
            tp = t

        fp.write('\end{document}')
        fp.close()            
        
        # Compile a pdf from the LaTex file.
        subprocess.call(['pdflatex',ftex],inline=True,shell=True)
        f.close()  


    def readegy(self,fegy=''):
        ''' 
        Read in the .egy file and convert the energies to MHz. Also read in the .cat
        file.
        '''
        if (len(fegy)<1):
            fegy=self.fegy
        self.fegy = fegy
        c = 299792458. 		# speed of light in m/s
        ceq = spfit.spfitPrediction(fname = fegy)
        ceq.read_egy()
        ceq.read_cat()
        self.egy_mhz = ceq.egy_content['egy']*c*100./1.e6
        self.egy_qn = ceq.egy_content['qn']
        self.cat_qn = ceq.cat_content['qn']
        self.cat_mhz = ceq.cat_content['freq']
        self.cat_logint = ceq.cat_content['lgint']

    def allowedtrans(self, J = 2, ka = 0, kc = 2):
        '''
        Looks for allowed transitions from the input. 

        Parameters
        ----------
        J (int): J quantum number of requested level (default = 2)
        ka (int): ka quantum number of requested level (default = 0)
        kc (int): kc quantum number of requested level (default = 2)
        
        
        Returns
        -------
        atype
        btype
        ctype
        
        Notes
        -----
        none
        
        '''
        if (len(self.egy_mhz)<3):
            self.readegy()
            
        Jto=J-1
        if (Jto<0):
            Jto=0
        
        atype = np.array([])
        btype = np.array([])
        ctype = np.array([])
        lvls = []
        Jto = Jto-1
        while (Jto<J+1):
            Jto = Jto + 1
            for kato in np.arange (0,Jto+1):
                kcto = Jto - kato
                if (kcto < Jto+1):
                    if (len(lvls) > 1):
                        lvls = np.vstack((lvls,[Jto,kato,kcto]))
                    if (len(lvls) < 1):
                        lvls = [Jto,kato,kcto]
                kcto = Jto + 1 - kato
                if (kcto < Jto+1):
                    if (len(lvls) > 1):
                        lvls = np.vstack((lvls,[Jto,kato,kcto]))
                    if (len(lvls) < 1):
                        lvls = [Jto,kato,kcto]   
        for x in lvls:
            Jt = x[0]
            kat = x[1]
            kct = x[2]
            for y in np.arange(len(self.egy_mhz)):
                Je = int(self.egy_qn[y][0:3])
                kae = int(self.egy_qn[y][3:7])
                kce = int(self.egy_qn[y][7:10])
                egt = 0.0
                if ((Je == Jt) and (kae == kat) and (kce == kct)):
                    egt = self.egy_mhz[y]
                    break
            if ((abs(Jt-J)<2) and (np.mod(abs(kat-ka),2)==0) and (np.mod(abs(kct-kc),2)!=0)):
                if (len(atype) > 1):
                    atype = np.vstack((atype,[Jt, kat, kct, egt]))
                if (len(atype) < 1):
                    atype = [Jt,kat,kct, egt]
            elif ((abs(Jt-J)<2) and (np.mod(abs(kat-ka),2)!=0) and (np.mod(abs(kct-kc),2)!=0)):
                if (len(btype) > 1):
                    btype = np.vstack((btype,[Jt, kat, kct, egt]))
                if (len(btype) < 1):
                    btype = [Jt,kat,kct,egt]            
            elif ((abs(Jt-J)<2) and (np.mod(abs(kat-ka),2)!=0) and (np.mod(abs(kct-kc),2)==0)):
                if (len(ctype) > 1):
                    ctype = np.vstack((ctype,[Jt, kat, kct, egt]))
                if (len(ctype) < 1):
                    ctype = [Jt,kat,kct,egt]
        return [atype,btype,ctype]

    def checktranstype(self,Jl,kal,kcl,Ju,kau,kcu):
        '''
        Determines whether a transition is a-, b-, or c-type or forbidden (f).
        
        return type = a, b, c, or f (string), where f means forbidden.
        '''
        if (abs(Jl-Ju)<2) and np.mod(abs(kal-kau),2)==0 and np.mod(abs(kcu-kcl),2)!=0:
            return 'a'
        elif (abs(Jl-Ju)<2) and np.mod(abs(kal-kau),2)!=0 and np.mod(abs(kcu-kcl),2)!=0:
            return 'b'
        elif (abs(Jl-Ju)<2) and np.mod(abs(kal-kau),2)!=0 and np.mod(abs(kcu-kcl),2)==0:
            return 'c'
        else:
            return 'f'

    def getlevelenergy(self,J,ka,kc):
        '''
        Search for and return the energy of a level as given in the .egy file.
        '''
        
        if (len(self.egy_mhz)<3):
            self.readegy()        

        for y in np.arange(len(self.egy_mhz)):
            Je = int(self.egy_qn[y][0:3])
            kae = int(self.egy_qn[y][3:7])
            kce = int(self.egy_qn[y][7:10])
            if (abs(Je-J)<.5) and (abs(kae-ka)<0.5) and (abs(kce-kc)<0.5):
                return self.egy_mhz[y]

    def searchintensity(self,quantumnumbers):
        '''
        Get the intensities for a transtion from the .cat file.
        
        quantnumbers is assumed to have the format:
        [[Jlower, kalower, kclower, Jupper, kaupper, kcupper, transition energy]]
        '''
        
        if (len(self.cat_qn)<3):
            self.readegy() 
        
        for x in quantumnumbers:
            Jl = int(x[0])
            kal = int(x[1])
            kcl = int(x[2])
            Ju = int(x[3])
            kau = int(x[4])
            kcu = int(x[5])
            for y in np.arange(len(self.cat_qn)):
                Ji = int(self.cat_qn[y][0:2])
                kai = int(self.cat_qn[y][2:4])
                kci = int(self.cat_qn[y][4:6])
                Jf = int(self.cat_qn[y][12:14])
                kaf = int(self.cat_qn[y][14:16])
                kcf = int(self.cat_qn[y][16:18])
                #print Jl,kal,kcl,Ju,kau,kcu
                #print Ji,kai,kci,Jf,kaf,kcf
                if Jl==Ji and kal==kai and kcl==kci and Ju==Jf and kau==kaf and kcu==kcf:
                    return self.cat_logint[y]
                elif Jl==Jf and kal==kaf and kcl==kcf and Ju==Ji and kau==kai and kcu==kci:
                    return self.cat_logint[y]
                elif (y==len(self.cat_qn)):
                    return -100.
                
            
    def searchlisdrv(self,drvt,tl,xl,xu,listenmax,listenmin,drivemax,drivemin):
        '''
        To Do
        - Add intensity information
        - Duplicate cycle checking and removing
        
        Checks if potential drive and listen transitions have the appropriate types
        and transition energies within the requested limits. It puts found cycles
        into the appropriate attribute.
        
        Please note: For three wave mixing cycle searches, it is expected that the
        transition type of the twist (xl -> xu) should differ from that of the lower 
        twist level to the inputted rotational levels in tl (xl -> tl), and that the
        drive type transitions (drvt) should be of the third transition type. For example,
        if (xu -> xl) is a-type and the drvt = 'b', then all (xl -> tl) should be c-type.
        

        Inputs
        ------
        drvt        : the desired type for the drive transition, a, b, or c (str)
        
        tl          : array containing quantum numbers and energy of a rotational level.
                      Note that tl can contain multiple rows.
                        tl = [[J, ka, kc, level energy]]
                        
        xl          : array containing the quantum numbers of the lower level of the
                      desired twist transition as well as the transition energy of
                      the twist transition. xl is a single row array.
                        xl = [transition energy, J, ka, kc]
                        
        xu          : array containing the quantum numbers of the upper level of the
                      desired twist transition. xu is a single row array.
                        xu  = [J, ka, kc]
                        
        listenmax   : The *max and *min variables are scalars that set the transition 
        listenmin     energy limits for searching for a cycle.
        drivemax
        drivemin
               
        Returns
        -------        
        None

        '''
        twst = np.array([])
        lstn = np.array([])
        drve = np.array([])

        for y in tl:
            egl = y[3]-self.getlevelenergy(xl[1],xl[2],xl[3])
            if abs(egl) <= listenmax and abs(egl) >= listenmin:
                drvtyp = self.checktranstype(y[0],y[1],y[2],xu[0],xu[1],xu[2])
                if drvtyp == drvt:
                    # Check if transition energy is within limits
                    egd = self.getlevelenergy(y[0],y[1],y[2])-self.getlevelenergy(xu[0],xu[1],xu[2])
                    if (abs(egd) <= drivemax) and abs(egd) >= drivemin:
                        # Cycle found, get the intensities.
                        inttw = self.searchintensity(quantumnumbers = [[xl[1],xl[2],xl[3],xu[0],xu[1],xu[2],xl[0]]])
                        intls = self.searchintensity(quantumnumbers = [[xl[1],xl[2],xl[3],y[0],y[1],y[2],egl]])
                        intdr = self.searchintensity(quantumnumbers = [[y[0],y[1],y[2],xu[0],xu[1],xu[2],egd]])
                        
                        if (len(twst)>1):
                            twst = np.vstack((twst,[xl[1],xl[2],xl[3],xu[0],xu[1],xu[2],xl[0],inttw]))
                            lstn = np.vstack((lstn,[xl[1],xl[2],xl[3],y[0],y[1],y[2],egl,intls]))
                            drve = np.vstack((drve,[y[0],y[1],y[2],xu[0],xu[1],xu[2],egd,intdr]))
                        if len(twst) < 2:
                            twst = [xl[1],xl[2],xl[3],xu[0],xu[1],xu[2],xl[0],inttw]
                            lstn = [xl[1],xl[2],xl[3],y[0],y[1],y[2],egl,intls]
                            drve = [y[0],y[1],y[2],xu[0],xu[1],xu[2],egd,intdr]
        # Put the found cycles, if any, into the appropriate attributes.
        if (len(self.twists)>1) and len(twst) > 2:
            self.twists = np.vstack((self.twists,twst))                    
            self.listens = np.vstack((self.listens,lstn))
            self.drives = np.vstack((self.drives,drve))
        elif (len(self.twists)<2) and len(twst) > 2:
            self.twists = twst
            self.listens = lstn
            self.drives = drve

    def threewavesearch(self,rfmin=50.,rfmax=560.,drivemin=2000.,drivemax=6000.,
                        listenmin=2000.,listenmax=6560.,Jmin=0,Jmax=7):
        '''
        Looks for threewave transitions. 

        Parameters
        ----------
        rfmin, rfmax         : Minimum and maximum RF (twist) frequencies
        drivemin, drive max  : Minimum and Maximum drive frequencies.
        listenmin, listenmax : Minimum and Maximum listen frequencies.
        Jmin, Jmax           : Minimum and Maximum J quantum numbers to consider in cycles.
        
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        # First, check if the .egy file has been read.
        if (len(self.egy_mhz)<3):
            self.readegy()
         
        self.rfmin = rfmin
        self.rfmax = rfmax
        self.dmin = drivemin
        self.dmax = drivemax
        self.lmin = listenmin
        self.lmax = listenmax
        self.Jmin = Jmin
        self.Jmax = Jmax
        
        # Now search for possible twist transitions.
        twists = np.array([])
        for x in np.arange(len(self.cat_mhz)):
             if ((self.cat_mhz[x] >= rfmin) and (self.cat_mhz[x] <= rfmax) 
             and (int(self.cat_qn[x][0:2]) >= Jmin) 
             and (int(self.cat_qn[x][12:14]) <= Jmax)):
                 tegy = self.cat_mhz[x]
                 tJi = int(self.cat_qn[x][0:2])
                 tkai = int(self.cat_qn[x][2:4])
                 tkci = int(self.cat_qn[x][4:6])
                 tJf = int(self.cat_qn[x][12:14])
                 tkaf = int(self.cat_qn[x][14:16])
                 tkcf = int(self.cat_qn[x][16:18])
                 if (len(twists)>1):
                     twists = np.vstack((twists,[tegy,tJi,tkai,tkci,tJf,tkaf,tkcf,self.cat_logint[x]]))
                 elif (len(twists)<1):
                     twists = [tegy,tJi,tkai,tkci,tJf,tkaf,tkcf,self.cat_logint[x]]
        self.cantwists = twists
        
        for x in twists:
            # check twist transition type
            trtyp = self.checktranstype(x[1],x[2],x[3],x[4],x[5],x[6])

            # Look for transitions from twist levels
            [atl,btl,ctl] = self.allowedtrans(x[1],x[2],x[3])
            [atu,btu,ctu] = self.allowedtrans(x[4],x[5],x[6])
            atl = np.vstack((atl,atu))
            btl = np.vstack((btl,btu))
            ctl = np.vstack((ctl,ctu))
            # For a-type twists, give either the b or c-type allowed transitions and look
            # for c or b-type drives, respectively.
            if trtyp == 'a':
                # For b-type listens
                self.searchlisdrv('c',btl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)
                # for c-type listens
                self.searchlisdrv('b',ctl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)

            # For b-type twists, give either the a or c-type allowed transitions and look
            # for c or a-type drives, respectively.
            if trtyp == 'b':
                # For a-type listens
                self.searchlisdrv('c',atl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)
                # for c-type listens
                self.searchlisdrv('a',ctl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)

            # For c-type twists, give either the b or a-type allowed transitions and look
            # for a or b-type drives, respectively.
            if trtyp == 'c':
                # For b-type listens
                self.searchlisdrv('a',btl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)
                # for c-type listens
                self.searchlisdrv('b',atl,x[0:4],x[4::1],listenmax,listenmin,drivemax,drivemin)

    def plot(self,J1=1,ka1=0,kc1=1,J2=2,ka2=0,kc2=2,J3=2,ka3=1,kc3=1,cycle_num=1e9):
        '''
        Plot the cycle specified by the quantum numbers or the cycle_num param.
        The quantum numbers are given as groups of 3 as Ji,kai,kci for the 
        i = 1, 2, and 3 levels. Alteratively, the parameter cycle_num may be given,
        which corresponds to the index of the desired cycle contained in self.twists,
        self.drives, and self.listens.
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        '''
        import pylab as pl

        if cycle_num > -1 and cycle_num < 1e9:
            if len(self.drives) < 2:
                print 'No cycles have been searched for and found yet. Please try a'
                print 'cycle search, then come back here to try plotting again.'
                print 'Alternatively, you can directly input the quantum numbers'
                print 'of the cycle you are interested in.'
                return
                
        if (cycle_num > -1 and cycle_num < len(self.drives)):
            d = self.drives[cycle_num]
            t = self.twists[cycle_num]
            
            [J1,ka1,kc1,J2,ka2,kc2] = d[0:6]
            
            if all(d[0:3]==t[0:3]):
                [J3,ka3,kc3]=t[3:6]
            else:
                [J3,ka3,kc3]=t[0:3]
        elif (cycle_num < -1) or (cycle_num < 1e9 and cycle_num > len(self.drives)):
            print 'Requested cycle_num = ',cycle_num,' is out of the range. Please'
            print 'try again.'
            return

        fntsz = 24        
        
        E1 = self.getlevelenergy(J1,ka1,kc1)
        E2 = self.getlevelenergy(J2,ka2,kc2)
        E3 = self.getlevelenergy(J3,ka3,kc3)
        
        Emin=min([E1,E2,E3])
        Emax=max([E1,E2,E3])
        
        [a1,b1,c1]=self.allowedtrans(J1,ka1,kc1)
        a1 = np.vstack((a1,b1))
        a1 = np.vstack((a1,c1))

        [a2,b2,c2]=self.allowedtrans(J2,ka2,kc2)
        a1 = np.vstack((a1,a2))
        a1 = np.vstack((a1,b2))
        a1 = np.vstack((a1,c2))

        [a3,b3,c3]=self.allowedtrans(J3,ka3,kc3)
        a1 = np.vstack((a1,a3))
        a1 = np.vstack((a1,b3))
        a1 = np.vstack((a1,c3))
        
        als = a1[a1[:,3].argsort()]

        ci = 0
        for x in als:
            if ci > 1:            
                if any(alnd[ci-1,:] != x):
                    if all(x[0:3]==[J1,ka1,kc1]) or all(x[0:3]==[J2,ka2,kc2]) or all(x[0:3]==[J3,ka3,kc3]):
                        continue
                    else:
                        alnd = np.vstack((alnd,x))
                        ci=ci+1
            if ci == 1:
                if all(x[0:3]==[J1,ka1,kc1]) or all(x[0:3]==[J2,ka2,kc2]) or all(x[0:3]==[J3,ka3,kc3]):
                    continue
                else:
                    if any(alnd != x):
                        alnd = np.vstack((alnd,x))
                        ci = ci+1
            if ci == 0:
                if all(x[0:3]==[J1,ka1,kc1]) or all(x[0:3]==[J2,ka2,kc2]) or all(x[0:3]==[J3,ka3,kc3]):
                    continue
                else:                
                    alnd = x
                    ci = 1
                
        [m2,m3] = [0.,0.]
        if (abs(E1-E2) < 500. or abs(E2-E3) < 500.):
            m2 = -1.
        if (abs(E1-E3) < 500.):
            m3 = -1.
        L1str = '$'+str(int(J1))+'_{'+str(int(ka1))+str(int(kc1))+'}$'
        L2str = '$'+str(int(J2))+'_{'+str(int(ka2))+str(int(kc2))+'}$'
        L3str = '$'+str(int(J3))+'_{'+str(int(ka3))+str(int(kc3))+'}$'

        self.figure = pl.figure(num=1,figsize=(15,20),dpi=300)
        pl.title('Three-wave-mixing cycle for '+L1str+', '+L2str+', '+L3str+' from '+self.fegy,fontsize=fntsz*.8)
        pl.plot([1.,2.],[E1,E1],color='black')
        pl.text(2.1,E1,L1str,fontsize=fntsz)
        pl.plot([1.,2.],[E2,E2],color='black')
        pl.text(2.1+m2*1.3,E2,L2str,fontsize=fntsz)
        pl.plot([1.,2.],[E3,E3],color='black')
        pl.text(2.1+m3*1.3,E3,L3str,fontsize=fntsz)    
        pl.text(1.0,Emin-2000.,L1str+'$\ \leftrightarrow\ $'+L2str+', '+'$'+'{0:8.3f}'.format(abs(E1-E2))+'\ MHz$',fontsize = fntsz)
        pl.text(1.0,Emin-2600.,L2str+'$\ \leftrightarrow\ $'+L3str+', '+'$'+'{0:8.3f}'.format(abs(E3-E2))+'\ MHz$',fontsize = fntsz)
        pl.text(1.0,Emin-3200.,L3str+'$\ \leftrightarrow\ $'+L1str+', '+'$'+'{0:8.3f}'.format(abs(E1-E3))+'\ MHz$',fontsize = fntsz)
        
        [Ep1,Ep2,Ep3]=[-1000.,-1000.,-1000.] 
        for x in alnd:
            if (abs(x[3]-Emin) < 6500. or abs(x[3]-Emax) < 6500.):
                pl.plot([-1.0,0],[x[3],x[3]],color='black')
                if (abs(Ep1-x[3]) < 500. and abs(Ep2-x[3]) < 500. and abs(Ep3-x[3]) < 500.):
                    xt = -1.4
                elif (abs(Ep1-x[3]) < 500. and abs(Ep2-x[3]) < 500. and abs(Ep3-x[3]) > 500.):
                    xt = 0.5
                elif (abs(Ep1-x[3]) < 500. and abs(Ep2-x[3]) > 500. and abs(Ep3-x[3]) > 500.):
                    xt = -1.2
                else:
                    xt = 0.1
                Ep3 = Ep2
                Ep2 = Ep1
                Ep1 = x[3]
                pl.text(xt,x[3],'$'+str(int(x[0]))+'_{'+str(int(x[1]))+str(int(x[2]))+'}$',fontsize=fntsz)
        
        pl.xlim(-1.5,2.5)
        pl.ylim(Emin-6500.,Emax+6500.)
        pl.tick_params(labelbottom='off',bottom='off',top='off')
        pl.yticks(fontsize = fntsz*.8,rotation=90)
        pl.ylabel('$Level\ Energy\ (MHz)$',fontsize=fntsz)
        
    def plotlevels(self, lowlim = 2000.,highlim = 8500.):
        '''Make a plot of the energy levels involving transitions within the
        specified limits.
        
        '''
        import pylab as pl
        lvlp = []
        legp = []
        lvlt = []
        legt = []
        
        leneg = self.cat_mhz
        itn = -1
        for x in self.cat_qn:
            itn = itn + 1
            [Ji,kai,kci,Jf,kaf,kcf] = [int(x[0:2]),int(x[2:4]),int(x[4:6]),int(x[12:14]),int(x[14:16]),int(x[16:18])]
            if leneg[itn] > lowlim and leneg[itn] < highlim:
                egt1 = self.getlevelenergy(Ji,kai,kci)
                egt2 = self.getlevelenergy(Jf,kaf,kcf)
                if len(legt) > 0.5:
                    tcnt1 = 1
                    tcnt2 = 1
                    for y in legt:
                        if abs(egt1-y) < 0.01:
                            tcnt1 = 0
                        if abs(egt2-y) < 0.01:
                            tcnt2 = 0
                    if tcnt1 > 0.5:
                        lvlt = np.append(lvlt,str(Ji)+'$_{'+str(kai)+str(kci)+'}$')                
                        legt = np.append(legt,self.getlevelenergy(Ji,kai,kci))
                    if tcnt2 > 0.5:
                        lvlt = np.append(lvlt,str(Jf)+'$_{'+str(kaf)+str(kcf)+'}$')
                        legt = np.append(legt,self.getlevelenergy(Jf,kaf,kcf))
                if len(legt) < 0.5:
                    lvlt = np.append(lvlt,str(Ji)+'$_{'+str(kai)+str(kci)+'}$')
                    legt = np.append(legt,self.getlevelenergy(Ji,kai,kci))
                    lvlt = np.append(lvlt,str(Jf)+'$_{'+str(kaf)+str(kcf)+'}$')
                    legt = np.append(legt,self.getlevelenergy(Jf,kaf,kcf))
#        for x in np.arange(len(lvlt)):
#            tcnt = 0
#            if x > 1:
#                for y in legp:
                    #print y
#                    if (abs(y-legt[x]) < 0.01):
#                        tcnt = 1
#                    if tcnt < 1:
#                        lvlp = np.append(lvlp,lvlt[x])
#                        legp = np.append(legp,legt[x])
#            if x < 1:
        lvlp = lvlt
        legp = legt
        offp = np.zeros(len(legp))
        for x in np.arange(len(lvlp)):        
            nx = []
            for y in np.arange(len(lvlp)):
                if (abs(legp[x]-legp[y])<2000. and y>x):
                    nx = np.append(nx,x)
            for y in np.arange(len(nx)):
                offp[nx[y]] = y/10.
                
        self.figurelevels = pl.figure(num=2,figsize=(15,20),dpi=300)
        pl.clf()
        for x in np.arange(len(lvlp)):
            pl.plot([0,1],[legp[x],legp[x]],'black')
            pl.text(1.01+offp[x],legp[x],lvlp[x],fontsize = 24)
        pl.xlim(-.5,1.5)
        pl.ylim(min(legp)-1000.,max(legp)+1000.)
        pl.ylabel('Energy (MHz)')             
        pl.savefig(self.fegy+'.eps',dpi=300,format='eps')
        print len(legp)