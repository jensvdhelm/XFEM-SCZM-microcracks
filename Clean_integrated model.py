# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 10:05:04 2023

@author: Jens van der Helm, Arshdeep Singh Brar
@owner: TU Delft
"""

"""

XFEM-SCZM-microcracks

 

Copyright 2024 <Arshdeep Singh Brar, Jens van der Helm>

 

Licensed under the Apache License, Version 2.0 (the "License");

you may not use this file except in compliance with the License.

You may obtain a copy of the License at

 

    http://www.apache.org/licenses/LICENSE-2.0

 

Unless required by applicable law or agreed to in writing, software

distributed under the License is distributed on an "AS IS" BASIS,

WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and

limitations under the License.

"""

"""

To run the script use : abaqus cae noGUI=Mesh_generator.py

"""
# Function for generating parallel XFEM domains
def PointInDomain(point, line1, line2, isLinePerpendicular=False):
    '''
    Parameters
    ----------
    point : TUPLE
        (x, y) value of the coordinate to check
    line1 : TUPLE
        (m, c) the value of m and c in the description
        y = mx + c
        if perpendicular:
            FLOAT
            c for the equation x = c 
    line2 : TUPLE, optional
        (m, c) the value of m and c in the description
        y = mx + c
        if perpendicular:
            FLOAT
            c for the equation x = c
    isLinePerpendicular : BOOL
        To check if the Line is 90 deg. The default is False.

    Returns
    -------
    

    '''
    if isLinePerpendicular:
        
        c1 = line1
        c2 = line2
        
        x, y = point
        
        val = x
        
        if val>min(c1,c2) and val<=max(c1,c2):
            return True
        
        else:
            return False
        
    
    else:
        
        m1, c1 = line1
        m2, c2 = line2
        
        x, y = point
        
        if m1 != m2:
            raise ValueError('Lines are not parallel')
        
        val = y - m1*x
        
        if val>min(c1,c2) and val<=max(c1,c2):
            
            return True
            
        else:
            
            return False

######################################################################
################### ABAQUS scripting starts here #####################
######################################################################

# Importing ABAQUS modules
from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import assembly
import step
import material
import load
import interaction

# python modules
import random
import math
import numpy as np

# To ensure abaqus viewports stays stable
session.viewports['Viewport: 1'].setValues(displayedObject=None)

# Creating an abaqus model
compositeModel = mdb.Model(name='Composite')

#Parameters for the mesh 
Nply = 23 #Number of plies
Lx = 10.  #Length of the ply in x dir (mm)
Ly = 10.  #Length of the ply in y dir (mm)
nx = 40   #Number of elements in x dir/ply 
ny = 40   #Number of elements in y dir/ply
nz = 1   #Number of elements in z dir/ply
# nz can be converted into a list

d = 0.5  # Distance between two enriched XFEM zones in a ply 

# Input for Knudsen dynamic viscosity
kb = 1.380649*10**(-23) #Boltzmann constant, m2 kg s-2 K-1,
M = 2.016  #gas molar mass of the molecule, kg/kmol 2.016 for hydrogen
Na=6.022*10**26 #avogadro 1/kmol
m=M/Na #weight of one molecule
R=8314 # J/kmol.K
dm= 297*10**(-12) #molecule diameter, m 

# Fibre diameter in micrometer, for T700 is 7 micrometer
char_length_microm = [7, 1, 10, 100, 1000]
char_length_microm2 = list(map(float, char_length_microm))
char_length_m = np.array(char_length_microm2)*10**-6 # diameter per fibre in micrometer
r = 0.5*char_length_m[0]  # radius of cross-sectional area of fibre
RR = 8.95e-6 # closest separation between two fibre centres, value backwards calculated for TC1225
R = 0.5*RR

ff = (math.pi)/(2*math.sqrt(3))*(r/R)**2  # fibre volume fraction, = 0.5549 for TC1225

# Elastic Material properties
Ef = 240e3 # MPa. Modulus of T700 fibre
Em = 4100 # MPa. Modulus of LM-PAEK 
E1 = Ef*ff + Em*(1-ff)
E2 = 1/((ff/Ef) + ((1-ff)/Em))
E3 = E2
nu12 = 0.33 #Value for TC1225 @ 138 - 218 K
nu13 = nu12 #Value for TC1225 @ 138 - 218 K
nu23 = 0.32 # NLR value
G12 = 4.3e3 #N/mm2 (MPa) Value for TC1225
G13 = G12 #N/mm2 (MPa) Value for TC1225
G23 = E2/(2*(1+nu23)) #N/mm2 (MPa)
S12 = 152 #MPa, IPSS
rho = 1587 #kg/m3. Already the value for LM-PAEK/T700 at 293 K

# CTE
alpha11 = 0 #negligible
alpha22 = 3.63e-5 # (m/m)/K. Already the value for LM-PAEK/T700, T<Tg
alpha33 = 3.69e-5 # (m/m)/K. Already the value for LM-PAEK/T700, T<Tg
alpha12 = 0 # can be considered orthotropic
alpha13 = 0 # can be considered orthotropic
alpha23 = 0 # can be considered orthotropic
CrysShrink = 0.081 # m/m. Already the value for LM-PAEK/T700, crystallization shrinkage

# Thermal conductivity 293 K. Already the value for LM-PAEK/T700
Tk11= 6.4 # W/m/K
Tk12= 0 # can be considered orthotropic
Tk22= 0.81 # W/m/K
Tk13= 0 # can be considered orthotropic
Tk23= 0 # can be considered orthotropic
Tk33= 0.74 # W/m/K 

# Thermal conductivity 473 K Already the value for LM-PAEK/T700
wTk11= 9.0 # W/m/K
wTk12= 0 # can be considered orthotropic
wTk22= 1.04 # W/m/K
wTk13= 0 # can be considered orthotropic
wTk23= 0 # can be considered orthotropic
wTk33= 0.91 # W/m/K 

Ft_mean = 86 #Mean value of transverse matrix strength (N/mm2)
m = 12 # Shape of Weibull curve taken from Grogan
beta = 0 #non-linearity of the shear stress–shear strain relation, which is zero for a linear behavior

GIc = 2.1 #Fracture energy under mode I (N/mm) Value for TC1225
GIIc = 2.6 #Fracture energy under mode II (N/mm) Value for TC1225
GIIIc = GIIc #Fracture energy under mode III (N/mm)

#Variables for different conditions 
p_constant_bar= [4, 500, 700] # in bar
T_var_Celcius= [-253, -220, 20] # in Celsius
p_constant_pascal= [x*10**5 for x in p_constant_bar]
T_var_Kelvin = [x+273.15 for x in T_var_Celcius]

#Thin ply lay-up
plyOrientations = [30., -30., 90., 30., -30., 90., 30., -30., 90., 30., -30., 90., -30., 30., 90., -30., 30., 90., -30., 30., 90., -30., 30., ] #List specifying ply Orientation
thickness = [0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110, 0.110] #mm, List specifying thickness of each ply

# Dynamic viscosity 
mu=[(2/(3*math.sqrt(math.pi)))*(math.sqrt(M*R*x)/(math.pi*Na*dm**2)) for x in T_var_Kelvin]

# Mean free path
labda = []
labda_nm = []
for k in range(len(mu)):
    labdazz=(mu[k]/p_constant_pascal[k])*math.sqrt(math.pi/2)*(math.sqrt(R*T_var_Kelvin[k]/M)) #if viscosity in Pa.s, pressure in Pa, R in (kgm2/s2)/kmolK, T in K,molecular mass in kg/kmol, then mean free path is in m
    labda_nmzz = labdazz*10**9 #nm
    labda.append(labdazz)
    labda_nm.append(labda_nmzz)

# Calculating the total thickness of the laminate
tot_thickness = 0
for t in thickness:
    tot_thickness = tot_thickness + t

# Calculating the shear strength of the matrix per ply
for i in range(Nply): 
    if thickness[i] <= 0.120:
        if beta != 0:
            if i == 0: # Thin, outer ply
                phi = (24*GIIc)/(math.pi*thickness[i])
            elif i == Nply-1: #Thin, outer ply
                phi = (24*GIIc)/(math.pi*thickness[i])
            else: #Thin, embedded ply
                phi = (48*GIIc)/(math.pi*thickness[i])
               
            Fs = math.sqrt(((1 + beta*phi*G12**2)**.5 - 1)/(3*beta*G12))    #In situ shear strength
            
        elif beta == 0:
            if i == 0: #Thin, outer ply
                Fs = 2*math.sqrt((G12*GIIc)/(math.pi*thickness[i]))
            elif i == Nply-1: #Thin, outer ply
                Fs = 2*math.sqrt((G12*GIIc)/(math.pi*thickness[i]))
            else: #Thin, embedded ply
                Fs = math.sqrt((8*G12*GIIc)/(math.pi*thickness[i])) #In situ shear strength
    if thickness[i] > 0.120:
        if beta != 0:
            if i == 0: # Thick, outer ply
                phi = (24*GIIc)/(math.pi*thickness[i])
            elif i == Nply-1: #Thick, outer ply
                phi = (24*GIIc)/(math.pi*thickness[i])
            else: #Thick, embedded ply
                phi = (12*S12**2)/(G12) + (72*beta*S12**4)/(4)
               
            Fs = math.sqrt(((1 + beta*phi*G12**2)**.5 - 1)/(3*beta*G12))    #In situ shear strength
            
        elif beta == 0:
            if i == 0: #Thick, outer ply
                Fs = 2*math.sqrt((G12*GIIc)/(math.pi*thickness[i]))
            elif i == Nply-1: #Thick, outer ply
                Fs = 2*math.sqrt((G12*GIIc)/(math.pi*thickness[i]))
            else: #Thick, embedded ply
                Fs = math.sqrt(2)*S12 #In situ shear strength

# Implementing permeability variables
e = 0.03  # Void ratio, JUST AN EXAMPLE YET! -> deduct it from micro-CT scan
leakrate = 1.07e-16 # m^2/s, taken from permeabilityresults thin ply LM-PAEK
k11= leakrate/tot_thickness # m/s
k12= leakrate/tot_thickness # m/s
k22= leakrate/tot_thickness # m/s
k13= leakrate/tot_thickness # m/s
k23= leakrate/tot_thickness # m/s
k33= leakrate/tot_thickness # m/s
gamma_w = 811.95 #N/m^3 specific weight of the wetting liquid = gaseous hydrogen @293K, 1 bara
# beta_e = #velocity coefficient, only needed for Forchheimers law

# Knudsen fitting parameters
Kn05 = 4.5
n = 5

ko = k11 #m^2, intrinsic permeability (actually diffusion) for a TC1225

# Defining the Knudsen number
Kn = []
f = []
dimensionless =[]
kapp = []
for j in range(len(char_length_microm)):
    for i in range(len(mu)):
        Knzz= labda[i]/char_length_m[j] #Knudsen number
        Kn.append(Knzz)
        fzz = 1/(1+(Knzz/Kn05)**(n-1)) #weighing coefficient
        f.append(fzz)
        dimensionlesszz = (1+4*Knzz)*fzz + (64*Knzz*(1-fzz))/(3*math.pi) #dimensionless permeability model 
        dimensionless.append(dimensionlesszz)
        kappzz = dimensionlesszz/ko # dimensionless apparent permeability
        kapp.append(kappzz)
        
        if Knzz <= 0.1:
            print("Diffusion by viscous/slip flow") 
        elif Knzz >= 10:
            print("Diffusion by Knudsen flow")
        else:
            print("Diffusion by transitional flow")
    

#Initiating Sketch of the frame in ABAQUS
plySketch = compositeModel.ConstrainedSketch(name='FrameSketch', sheetSize=30)

# Creating a rectangle 
plySketch.rectangle(point1=(-Lx/2, -Ly/2), point2=(Lx/2, Ly/2))

# Creating part and meshing the region  
for i in range(Nply):
    
    #Naming the part
    partname = 'Ply_' + str(i+1)
    #Initiating the part
    thispart = compositeModel.Part(name=partname, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    #Extruding the part
    thispart.BaseSolidExtrude(sketch=plySketch, depth=thickness[i])
    
    #Storing edges and cell in a local variable
    Edges = thispart.getFeatureEdges(name='Solid extrude-1')
    Cells = thispart.getFeatureCells(name='Solid extrude-1')
    
    #Create set for the whole ply
    plyset = thispart.Set(name='ply-'+str(i+1), cells=Cells)
    
    matOrientation = thispart.MaterialOrientation(orientationType=SYSTEM, region=plyset, 
                                                  axis=AXIS_3, angle=plyOrientations[i],
                                                  stackDirection=STACK_3)
    # Surface sets for delamination zones
    # For the first ply and the last ply only surface to be used is the top and bottom surface respectively
    if i == 0:
        surfName = 'Ply-' + str(i+1) + '-topSurf'
        thispart.Surface(name=surfName, 
                         side1Faces=thispart.faces.findAt(((0., 0., thickness[i]),)))
        
    elif i == Nply-1:
        surfName = 'Ply-' + str(i+1) + '-bottomSurf'
        thispart.Surface(name=surfName,
                         side1Faces=thispart.faces.findAt(((0., 0., 0.),)))
        thispart.Set(name='xySym', faces=(thispart.faces.findAt(((0., 0., thickness[i]),)),))
        
    else:
        surfName1 = 'Ply-' + str(i+1) + '-topSurf'
        surfName2 = 'Ply-' + str(i+1) + '-bottomSurf'
        thispart.Surface(name=surfName1, 
                         side1Faces=thispart.faces.findAt(((0., 0., thickness[i]),)))
        thispart.Surface(name=surfName2,
                         side1Faces=thispart.faces.findAt(((0., 0., 0.),)))
    
    #set for faces to assign boundary conditions
    thispart.Set(name='ply-'+str(i+1)+'-xminus', 
                 faces=(thispart.faces.findAt(((-Lx/2, 0., thickness[i]/2),)),))
    
    thispart.Set(name='ply-'+str(i+1)+'-xplus', 
                 faces=(thispart.faces.findAt(((Lx/2, 0., thickness[i]/2),)),))
    
    thispart.Set(name='ply-'+str(i+1)+'-yminus', 
                 faces=(thispart.faces.findAt(((0., -Ly/2, thickness[i]/2),)),))
    
    thispart.Set(name='ply-'+str(i+1)+'-yplus', 
                 faces=(thispart.faces.findAt(((0., Ly/2, thickness[i]/2),)),))
    
    # Selecting the region to be meshed
    thisregion = regionToolset.Region(cells=Cells)
    
    #Assigning the element type for the mesh
    thispart.setElementType(regions=thisregion, elemTypes=(mesh.ElemType(elemCode=DC3D8),))
    
    #Seeding edges according to the number of elements 
    thispart.seedEdgeByNumber(edges=(thispart.edges.findAt((0, Ly/2., 0.)),
                                     thispart.edges.findAt((0, Ly/2., thickness[i])),
                                     thispart.edges.findAt((0, -Ly/2., 0.)),
                                     thispart.edges.findAt((0, -Ly/2., thickness[i]))),
                              number=nx)
    thispart.seedEdgeByNumber(edges=(thispart.edges.findAt((Lx/2., 0., 0.)),
                                     thispart.edges.findAt((Lx/2., 0., thickness[i])),
                                     thispart.edges.findAt((-Lx/2., 0., 0.)), 
                                     thispart.edges.findAt((-Lx/2., 0., thickness[i]))),
                              number=ny)
    thispart.seedEdgeByNumber(edges=(thispart.edges.findAt((Lx/2., Ly/2., thickness[i]/2)), 
                                     thispart.edges.findAt((Lx/2., -Ly/2., thickness[i]/2)),
                                     thispart.edges.findAt((-Lx/2., Ly/2., thickness[i]/2)),
                                     thispart.edges.findAt((-Lx/2., -Ly/2., thickness[i]/2))),
                              number=nz)
    thispart.generateMesh()
    
for i in range(Nply):
    #Note: uncapitatlized plural words indicate already existing data 
    # Storing the part object in thispart
    thispart = compositeModel.parts['Ply_'+str(i+1)]
    
    # This loop assigns material to each element
    for element in thispart.elements:
        # Creating element sequence datatype to be later used for Region command
        MeshElem = thispart.elements.sequenceFromLabels(labels=(element.label,))
        # Creating a region out of current element
        elemRegion = regionToolset.Region(elements=MeshElem)
        
        # Creating name for the material and section
        MatName = 'TC1225-' + str(i+1) + '-' + str(element.label) 
        SectName = 'Section_CFRP-' + str(i+1) + '-' + str(element.label)
        
        # Creating a material
        compositeModel.Material(name=MatName)
        
        # Generating random strength distribution
        Ft = random.weibullvariate(alpha=Ft_mean, beta=m)
        
        
#Implementation of randomly distributed voids (ellipsoids)
    # mux =  # mean, from CT-scan
    # stdx = # standard deviation, from CT-scan
    # muy = 
    # stdy = 
    # muz = 
    # stdz = 
    
    # NrofV = round((3*Vv*Lx*Ly*tot_thickness)/(4*math.pi*0.5*mux*0.5*muy*0.5*muz))  # number of voids, making use of ellipse volume formula
    
    # from scipy.stats import norm  # generate random numbers from N(0,1). Normal distribution / Gaussian distribution
    # data_normalx = norm.rvs(size=NrofV,loc=mux,scale=stdx)  # 'size is the number of random variates'
    # data_normaly = norm.rvs(size=NrofV,loc=muy,scale=stdy) 
    # data_normalz = norm.rvs(size=NrofV,loc=muz,scale=stdz)
    # void_sizes = zip(data_normalx,data_normaly,data_normalz) #  set of x-, y-, and z-size per void
    # print(void_sizes)

    #Alternative void implementation by spherical voids
     #   compositeModel.rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
      #      mdb.models['Model-1'].rootAssembly.instances['sphere-1'],
        
        
        # Temperature dependent conductivity    
        compositeModel.materials[MatName].Conductivity(table=((Tk11, Tk12, Tk22, Tk13, Tk23, Tk33, 293.0), 
                                                              (wTk11, wTk12, wTk22, wTk13, wTk23, wTk33, 473.0)),
                                                                   temperatureDependency=ON, type=ANISOTROPIC)
    
        compositeModel.materials[MatName].Density(table=((rho, ), ))
        
        # Temperature independent CTE (anisotropic) (only temperature independent below Tg = 313K)
        compositeModel.materials[MatName].Expansion(table=((alpha11, alpha22, alpha33, alpha12, alpha13, 
                                                            alpha23), ), type=ANISOTROPIC)
        
        compositeModel.materials[MatName].Elastic(type=ENGINEERING_CONSTANTS, 
                                                  table=((E1, E2, E3, nu12, nu13, nu23, 
                                                         G12, G13, G23),))
        
        # Specific heat for TC1225 below Tg
        compositeModel.materials[MatName].SpecificHeat(table=((820.0, 20.0), 
                                                                          (982.0, 50.0), (1144.0, 100.0), (1306.0, 150.0), (1468.0, 200.0), (1500.0, 
                                                                                                                                           230.0)), temperatureDependency=ON)
        
        # Implementing permeability by Darcy's law
#        if Kn[1,1] <= 0.001:
        compositeModel.materials[MatName].Permeability(type=ANISOTROPIC, specificWeight=gamma_w, 
                                               temperatureDependency=ON, table=((k11, k12, 
                                                k22, k13, k23, k33, e, T_var_Kelvin),))
    
        # Implementing permeability by Forchheimer's law
#        if Kn[1,1] >= 0.001:
#        compositeModel.materials[MatName].Permeability(type=ANISOTROPIC, specificWeight=gamma_w, 
#                                             temperatureDependency=ON, table=((k11, k12, 
#                                             k22, k13, k23, k33, e, T_var_Kelvin),))
#        compositeModel.material[MatName].permeability.VelocityDependence(table=((beta_e, e),))
#   
#        compositeModel.materials[MatName].MaxsDamageInitiation(direction=TMORI, 
#                                           table=((Ft, Fs, Fs),), tolerance=0.07)
        
        compositeModel.materials[MatName].maxsDamageInitiation.DamageEvolution(
                    type=ENERGY, softening=LINEAR, mixedModeBehavior=BK, power=1.45,
                    table=((GIc, GIIc, GIIIc),))
        
        compositeModel.materials[MatName].maxsDamageInitiation.DamageStabilizationCohesive(
                    cohesiveCoeff=1e-4)
        
        compositeModel.HomogeneousSolidSection(name=SectName, material=MatName, thickness=None)
        
        thispart.SectionAssignment(region=elemRegion, sectionName=SectName, offset=0., 
                                   offsetType=MIDDLE_SURFACE, offsetField='', 
                                   thicknessAssignment=FROM_SECTION)
        

# Setting up assembly for the model        
assembly = compositeModel.rootAssembly
#Assigning coordinate system to the assembly, cartesian in this case
assembly.DatumCsysByDefault(CARTESIAN)

#Assemblying one ply over the other
dist = -tot_thickness
for i in range(Nply):
    
    partName = 'Ply_' + str(i+1)
    
    thispart = compositeModel.parts[partName]
    
    thisInstance = assembly.Instance(name=partName+'-1', part=thispart, dependent=ON)
    
    thisInstance.translate(vector=(0.,0.,dist))
    
    dist = dist + thickness[i]


# Creating coupled temperature displacement step (thermal+mechanical load)    
#is performed when the mechanical and thermal solutions affect each other strongly and, therefore, must be obtained simultaneously
loadStep = compositeModel.CoupledTempDisplacementStep(deltmx=100.0, description=
    'Load_Step', initialInc=0.001, maxInc=0.001, maxNumInc=100000, minInc=1e-12
    , name='loadStep', nlgeom=ON, previous='Initial')

# Adding control parameters on the step
loadStep.control.setValues(allowPropagation=OFF, timeIncrementation=DEFAULT)

# Solution control for convergence
loadStep.control.setValues(discontinuous=ON)

# Creating delamination contact property
Interactionprop = compositeModel.ContactProperty('Delamination_prop')

Interactionprop.CohesiveBehavior(defaultPenalties=OFF,
                                 table=((1e8, 1e8, 1e8),))

Interactionprop.Damage(evolTable=((0.19, 0.79, 0.79),), evolutionType=ENERGY,
                       exponent=1.45, initTable=((51.7, 40.0, 40.0),), mixedModeType=BK,
                       useEvolution=ON, useMixedMode=ON, useStabilization=ON, 
                       viscosityCoef=1e-5)

# Creating contact interaction between plies for delamination
for i in range(Nply-1):
    
    instance1name = 'Ply_' + str(i+1) + '-1'
    surf1name = 'Ply-' + str(i+1) + '-topSurf'
    
    instance2name = 'Ply_' + str(i+2) + '-1'
    surf2name = 'Ply-' + str(i+2) + '-bottomSurf'
    
    compositeModel.SurfaceToSurfaceContactStd(adjustMethod=NONE, clearanceRegion=None,
                       createStepName='Initial', datumAxis=None,
                       initialClearance=OMIT, interactionProperty='Delamination_prop',
                       master=assembly.instances[instance1name].surfaces[surf1name], 
                       name='Delamination_layer'+str(i+1),
                       slave=assembly.instances[instance2name].surfaces[surf2name],
                       sliding=SMALL, thickness=ON)


# Set that can be used for boundary conditions or loading
FixedX = []
FixedY = []
LoadX = []
LoadY = []

for i in range(Nply):
    
    instanceName = 'Ply_' + str(i+1) + '-1'
    xfixedSet = 'ply-' + str(i+1) + '-xminus'
    yfixedSet = 'ply-' + str(i+1) + '-yminus'
    xloadSet = 'ply-' + str(i+1) + '-xplus'
    yloadSet = 'ply-' + str(i+1) + '-yplus'
    
    # Appending load and fixed sets of each ply
    FixedX.append(assembly.allInstances[instanceName].sets[xfixedSet])
    FixedY.append(assembly.allInstances[instanceName].sets[yfixedSet])
    LoadX.append(assembly.allInstances[instanceName].sets[xloadSet])
    LoadY.append(assembly.allInstances[instanceName].sets[yloadSet])
    
    ###### Here XFEM region assignments starts ########
    
    # Checks if the orientation is 90
    if plyOrientations[i] == 90.:
        isLinePerpendicular = True
        Nregions = int(Lx/d)
    else:
        isLinePerpendicular = False
        Nregions = int(Ly/d)
    
    # Initiating lines array 
    # Lines are used to define the array
    lines = []
        
    for j in range(Nregions + 1):
        
        if isLinePerpendicular:
            lines.append(-Lx/2 + j*d)
            
        else:
            lines.append((math.tan(plyOrientations[i]*math.pi/180), 
                          -Lx/2 + j*d/math.cos(plyOrientations[i]*math.pi/180)))
    XFEMRegions = {}
    
    for j in range(Nregions):
        XFEMRegions['Region'+str(j+1)] = []
    
    XFEMRegions['Region'+str(Nregions+1)] = []
        
    for element in assembly.allInstances[instanceName].elements:
        
        elemNodes = element.getNodes()
        
        xc = 0
        yc = 0
        
        #Calculating centroid of the element
        for node in elemNodes:
            xc = xc + node.coordinates[0]/8.
            yc = yc + node.coordinates[1]/8.
        
        elemCentroid = (xc, yc)
        
        elementFound = False
        
        for region in range(1, len(lines)):
            
            if PointInDomain(elemCentroid, lines[region-1], lines[region], isLinePerpendicular):
                
                XFEMRegions['Region'+str(region)].append(element.label)
                elementFound = True
                
                continue
        
        if not elementFound:
            XFEMRegions['Region' + str(Nregions+1)].append(element.label)
        
    print(XFEMRegions)
    
    for key, elems in XFEMRegions.items():
        
        if len(elems) != 0:
    
            assembly.Set(name=key+'-'+str(i+1), 
                    elements=assembly.allInstances[instanceName].elements.sequenceFromLabels(labels=elems))        
            XFEMregion = regionToolset.Region(
                elements=assembly.allInstances[instanceName].elements.sequenceFromLabels(labels= elems))
            assembly.engineeringFeatures.XFEMCrack(crackDomain=XFEMregion, name='Crack-'+str(i+1)+key)
    
    ###### XFEM region assignments ends ########


assembly.SetByBoolean(name='FixedX', sets=FixedX)
assembly.SetByBoolean(name='FixedY', sets=FixedY)
assembly.SetByBoolean(name='LoadX', sets=LoadX)
assembly.SetByBoolean(name='LoadY', sets=LoadY)

#XFixed boundary condition for FixedX
#YFixed boundary condition for FixedY
#SET implies constrained dof
compositeModel.DisplacementBC(name='XFixed', createStepName='Initial', region=
                              assembly.allSets['FixedX'], u1=SET, ur2=SET, ur3=SET)

compositeModel.DisplacementBC(name='YFixed', createStepName='Initial', region=
                              assembly.allSets['FixedY'], u2=SET, ur1=SET, ur3=SET)

compositeModel.DisplacementBC(name='XYSym', createStepName='Initial', region=
                              assembly.instances['Ply_'+str(Nply)+'-1'].sets['xySym'], u3=SET, ur1=SET, ur2=SET)

compositeModel.DisplacementBC(name='XLoad', createStepName='Load_Step', region=
                              assembly.allSets['LoadX'], u1=0.1)

compositeModel.DisplacementBC(name='YLoad', createStepName='Load_Step', region=
                              assembly.allSets['LoadY'], u2=0.1)


c1 = assembly.instances['Ply_1-1'].cells
c2 = assembly.instances['Ply_2-1'].cells
c3 = assembly.instances['Ply_3-1'].cells
c4 = assembly.instances['Ply_4-1'].cells
c5 = assembly.instances['Ply_5-1'].cells
c6 = assembly.instances['Ply_6-1'].cells
c7 = assembly.instances['Ply_7-1'].cells
c8 = assembly.instances['Ply_8-1'].cells
c9 = assembly.instances['Ply_9-1'].cells
c10 = assembly.instances['Ply_10-1'].cells
c11 = assembly.instances['Ply_11-1'].cells
c12 = assembly.instances['Ply_12-1'].cells
c13 = assembly.instances['Ply_13-1'].cells
c14 = assembly.instances['Ply_14-1'].cells
c15 = assembly.instances['Ply_15-1'].cells
c16 = assembly.instances['Ply_16-1'].cells
c17 = assembly.instances['Ply_17-1'].cells
c18 = assembly.instances['Ply_18-1'].cells
c19 = assembly.instances['Ply_19-1'].cells
c20 = assembly.instances['Ply_20-1'].cells
c21 = assembly.instances['Ply_21-1'].cells
c22 = assembly.instances['Ply_22-1'].cells
c23 = assembly.instances['Ply_23-1'].cells
assembly.Set(cells=c1[0:1]+c2[0:1]+c3[0:1]+c4[0:1]+c5[0:1]+c6[0:1]+c7[0:1]+c8[0:1]+c9[0:1]+c10[0:1]+c11[0:1]+c12[0:1]+c13[0:1]+c14[0:1]+c15[0:1]+c16[0:1]+c17[0:1]+c18[0:1]+c19[0:1]+c20[0:1]+c21[0:1]+c22[0:1]+c23[0:1], name='ALL')


# Implementing a predefined field for the desired temperature 
compositeModel.Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(20.0, ), name='Temperature_Predefined Field-1', region=
    assembly.allSets['ALL'])
  
# Implementing heat flux  
region = assembly.instances['Ply_1-1'].surfaces['Ply-1-bottomSurf']
compositeModel.SurfaceHeatFlux(createStepName='loadStep', distributionType=UNIFORM, magnitude=20.0, name='heatflux', region=region)
    
# Implementing convection boundary conditions    
compositeModel.FilmCondition(createStepName='loadStep', definition=
    EMBEDDED_COEFF, filmCoeff=0.001, filmCoeffAmplitude='', name='Int-1', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature=-273.0, surface=assembly.instances['Ply_1-1'].surfaces['Ply-1-bottomSurf']) 
    


Job = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Composite', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='TC1225', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

Job.writeInput()


            
    
                    
            
    
    
    





    
    


