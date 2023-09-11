import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.path as pltPath
import matplotlib.patches as pltPatches
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.coordinates import SkyOffsetFrame, ICRS
from astropy.io import ascii
from tabulate import *
'''
This program recieves a file with a list of RA_DEC in degrees or hms (simbad format) or galactic LON_LAT in degrees.
All in the same format. 
Returns 20 files 'FXX_extra.dat' with the RA_DEC visible in each 'XX' field.

::Examples::
160.9894113830400 -59.5476096936900
10 43 57.4587319296 -59 32 51.3948
287.41 -0.573805

Is necessary to create the "Output" folder at the location of the .py
Linux users will perhaps need to change the path ('\' to '/') 
'''

# G. Holgado 2023
########################################## INPUT
Archive='ListaCoordsPrueba_hms.txt'             # input file with coordinates; RA_DEC in degrees or hmsdms (simbad format i.e. spaces) or LON_LAT.
frameT='icrs'                                   # frame of coords: galactic or icrs ; ra&dec in deg is 'icrs'; ra&dec in hmsdms is 'icrs'; gal in deg is 'galactic'
extra='Proof'                                   # tag for files
ProjPlot='gal'                                  # Projection for plot (just for visual inspection): gal Aitoff
#########################################################################################################################################
####################################             NOT NECESSARY TO CONTINUE              #################################################
#########################################################################################################################################


########################################## VERTICES
def PintCuadrSkyCoord(ThisCenter,TypeElecc):
    center = ThisCenter
    Rango=1.45
    Mitad=Rango/2.
    Mitad_Rad=math.radians(Rango)/2.
    ##
    V1 = SkyCoord((ThisCenter.ra.degree-(Mitad/np.cos(ThisCenter.dec.rad-Mitad_Rad)))*u.deg,(ThisCenter.dec.degree-Mitad)*u.deg, frame='icrs')
    V2 = SkyCoord((ThisCenter.ra.degree+(Mitad/np.cos(ThisCenter.dec.rad-Mitad_Rad)))*u.deg,(ThisCenter.dec.degree-Mitad)*u.deg, frame='icrs')
    V3 = SkyCoord((ThisCenter.ra.degree+(Mitad/np.cos(ThisCenter.dec.rad+Mitad_Rad)))*u.deg,(ThisCenter.dec.degree+Mitad)*u.deg, frame='icrs')
    V4 = SkyCoord((ThisCenter.ra.degree-(Mitad/np.cos(ThisCenter.dec.rad+Mitad_Rad)))*u.deg,(ThisCenter.dec.degree+Mitad)*u.deg, frame='icrs')
    ##
    if TypeElecc=='gal':
        return (V1.galactic,V2.galactic,V3.galactic,V4.galactic)
    if TypeElecc=='Aitoff':
        return (V1.galactic.spherical,V2.galactic.spherical,V3.galactic.spherical,V4.galactic.spherical)
    


########################################## MUDEHaR POINTINGS DATA
A = SkyCoord(ra=[000.70, 013.21, 035.03, 036.11, 038.57, 043.39, 045.76, 056.74, 079.44, 083.85, 083.85, 085.17, 097.98, 100.22, 277.86, 295.73, 307.78, 308.29, 308.81, 313.70]*u.degree, dec=[+67.15, +56.63, +57.20, +62.05, +61.57, +60.46, +60.30, +24.07, +33.32, -04.60, -05.80, -01.97, +04.94, +09.81, -02.08, +23.49, +41.14, +40.74, +46.84, +43.81]*u.degree, frame='icrs')
A_gal = A.galactic

Coincid=[] # Final Vector with size 20. index from Archive visible in each pointing

DiccIApunt={'0':'000.70  +67.15  Berkeley_59',
'1':'013.21  +56.63  Pacman_Nebula',
'2':'035.03  +57.20  h_and_chi_Per',
'3':'036.11  +62.05  BD_+61_411_A_(O6.5_((f))z)',
'4':'038.57  +61.57  IC_1805_Heart',
'5':'043.39  +60.46  IC_1848_Soul',
'6':'045.76  +60.30  HD_18_326_S_(O6.5_V((f))z)',
'7':'056.74  +24.07  Pleiades',
'8':'079.44  +33.32  NGC_1893_C',
'9':'083.85  -04.60  42_Ori_Aa_(B0.7_V_+_B1:_V(n))',
'10':'083.85  -05.80  theta1_Ori_CaCb_(O7_Vp)',
'11':'085.17  -01.97  sigma_Ori_AaAb_C',
'12':'097.98  +04.94  NGC_2244_C',
'13':'100.22  +09.81  NGC_2264',
'14':'277.86  -02.08  W_40',
'15':'295.73  +23.49  NGC_6823',
'16':'307.78  +41.14  Cyg_OB2_CS',
'17':'308.29  +40.74  Cyg_OB2_CN',
'18':'308.81  +46.84  Berkeley_90_C',
'19':'313.70  +43.81  North_America_Nebula_CW'}

########################################## RA_DEC to try visible
Lista_Coords_Prueba=ascii.read(Archive,guess=True)
if len(Lista_Coords_Prueba[0])==2:
    Lista_Coords_Prueba_2RADEC=np.genfromtxt(Archive,dtype=str,delimiter=',')
    P = SkyCoord(Lista_Coords_Prueba_2RADEC,unit=(u.deg, u.deg), frame=frameT)
    FMT_Tab=(".11f",".11f")
if len(Lista_Coords_Prueba[0])!=2:
    Lista_Coords_Prueba_2hms=np.genfromtxt(Archive,dtype=str,delimiter=',')
    P = SkyCoord(Lista_Coords_Prueba_2hms, unit=(u.hourangle, u.deg), frame=frameT)
    FMT_Tab=(".0f",".0f",".11f",".0f",".0f",".11f")

P_gal = P.galactic

##################### GAL
if ProjPlot=='gal':
    fig, ax = plt.subplots(figsize=(6.5, 5.2),constrained_layout=True)
    cs = ax.scatter(A_gal.l.degree,A_gal.b.degree, s=34,zorder=1)
    ax.set_xlabel('Galactic longitude, $l$ [deg]')
    ax.set_ylabel('Galactic latitude, $b$ [deg]')
    for i,e in enumerate(A):
        temp_iP=[]
        ###PATCH OF FOV
        Recuadro = PintCuadrSkyCoord(e,ProjPlot)
        TheOnePath = pltPath.Path([[Recuadro[0].l.deg,Recuadro[0].b.deg],[Recuadro[1].l.deg,Recuadro[1].b.deg],
                                   [Recuadro[2].l.deg,Recuadro[2].b.deg],[Recuadro[3].l.deg,Recuadro[3].b.deg],
                                   [Recuadro[0].l.deg,Recuadro[0].b.deg]])
        TheOnePatch = pltPatches.PathPatch(TheOnePath,facecolor='gold',alpha=0.5) #PATCH TO PLOT 
        ax.add_patch(TheOnePatch)         
        ###Check_Distance 0.725 in RA and DEC
        center = e
        center.transform_to(SkyOffsetFrame(origin=center))
        Dist=P.transform_to(SkyOffsetFrame(origin=center))
        temp_iP.append(np.where(np.logical_and(abs(Dist.lat.deg)<0.725, abs(Dist.lon.deg)<0.725))[0])
        cs = ax.scatter(P_gal[temp_iP].l.deg,P_gal[temp_iP].b.deg, s=14,c='red',zorder=3,alpha=0.6)
        plt.annotate('F'+str(i+1).zfill(2),[e.galactic.l.deg-0.1,e.galactic.b.deg]) # 
        Coincid.append(temp_iP)
        Lista=True
        if Lista==True:
            print('F'+str(i+1).zfill(2)+': '+DiccIApunt[str(i)]               )
            print('-----------------------------------------------')
            print(Lista_Coords_Prueba[tuple(temp_iP)])            
            print('===============================================================')
            if len(temp_iP[0])!= 0:
                f = open('.\Output\F'+str(i+1).zfill(2)+'_'+extra+'.dat', 'w')                      
                f.write(tabulate(Lista_Coords_Prueba[tuple(temp_iP)],numalign='right',stralign='left',tablefmt='plain',floatfmt=FMT_Tab))
                f.close()
    plt.xlim(220,20)    
    #plt.ylim(3.75,5.75)
    plt.show() 

##################### Aitoff
if ProjPlot=='Aitoff':
    fig, ax = plt.subplots(figsize=(6.5, 5.2),constrained_layout=True,subplot_kw=dict(projection="aitoff")) ##aitoff, mollweide, hammer, lambert
    sph = A_gal.spherical
    cs = ax.scatter(-sph.lon.wrap_at(180*u.deg).radian,sph.lat.radian)
    ax.set_xlabel('Galactic longitude, $l$ [deg]')
    ax.set_ylabel('Galactic latitude, $b$ [deg]')
    for i,e in enumerate(A):
        temp_iP=[]
        ###PATCH OF FOV        
        Recuadro = PintCuadrSkyCoord(e,ProjPlot)
        TheOnePath = pltPath.Path([[-Recuadro[0].lon.wrap_at(180*u.deg).radian,Recuadro[0].lat.radian],[-Recuadro[1].lon.wrap_at(180*u.deg).radian,Recuadro[1].lat.radian],
                                   [-Recuadro[2].lon.wrap_at(180*u.deg).radian,Recuadro[2].lat.radian],[-Recuadro[3].lon.wrap_at(180*u.deg).radian,Recuadro[3].lat.radian],
                                   [-Recuadro[0].lon.wrap_at(180*u.deg).radian,Recuadro[0].lat.radian]])
        TheOnePatch = pltPatches.PathPatch(TheOnePath,facecolor='gold',alpha=0.5) #PATCH TO PLOT 
        ax.add_patch(TheOnePatch)         
        ###Check_Distance 0.725 in RA and DEC
        center = e
        center.transform_to(SkyOffsetFrame(origin=center))
        Dist=P.transform_to(SkyOffsetFrame(origin=center))
        temp_iP.append(np.where(np.logical_and(abs(Dist.lat.deg)<0.725, abs(Dist.lon.deg)<0.725))[0])
        if len(temp_iP)!= 0:
            sphP = P_gal.spherical
            cs = ax.scatter(-sphP[temp_iP].lon.wrap_at(180*u.deg).radian,sphP[temp_iP].lat.radian, s=14,c='red',zorder=2,alpha=0.6)#, vmin=1.5, vmax=2.5, cmap='twilight')                                
        Coincid.append(temp_iP)
        Lista=True
        if Lista==True:
            print('F'+str(i+1).zfill(2)+': '+DiccIApunt[str(i)]               )
            print('-----------------------------------------------')
            print(Lista_Coords_Prueba[tuple(temp_iP)])            
            print('===============================================================')
    ax.grid()
    plt.show() 

##########################################################
    


    
    
