print("*********  A python script to perform analaysis on an Enzo dataset using yt  *********")
print(" ")
print("Usage: Run the script as python yt_analysis.py FirstDataDump SecondDataDump PlotType (1 = slice, 2 = projection, 3 = phase, 4 = radial)")
print("       Run the script as python yt_analysis.py --h for more help")
#print("       Run the script as python yt_analysis.py --v for increased verbosity ")
#print("       The same as above but any value after sets the value of the Hubble constant, Omega Matter and Omega Vacuum, in that order")
print(" ")


#Packages
import yt
import math
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from yt.units import mp
from yt.units import kboltz
from numpy import linalg as LA
with np.errstate(divide='ignore'):
    np.float64(1.0) / 0.0
from yt.analysis_modules.level_sets.api import *
import matplotlib.cm as cm
from yt.utilities.physical_constants import gravitational_constant_cgs as G
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.data_objects.particle_filters import add_particle_filter

#Constants
pi = np.pi
HydrogenFractionByMass   = 0.76
DeuteriumToHydrogenRatio = 3.4e-5 * 2.0
HeliumToHydrogenRatio    = (1.0 - HydrogenFractionByMass) / HydrogenFractionByMass
SolarMetalFractionByMass = 0.01295

#Set sys argument to numberical argument 
argvs = sys.argv


#Data set
number_min = int(argvs[1])
number_max = int(argvs[2])

#Plot type
plot_type = ''
if int(sys.argv[3]) == 1:
    plot_type = 'Slice'
elif int(sys.argv[3]) == 2:
    plot_type = 'Proj'
elif int(sys.argv[3]) == 3:
    plot_type = 'Phase'
elif int(sys.argv[3]) == 4:
    plot_type = 'Radial'
else:
    print('Invalid')
    exit


#print(plot_switch)
#print(argvs[3])
#Slice = int(argvs[3])
#Proj = int(argvs[4])
#Phase = int(argvs[5])
#Radial = int(argvs[6])


#Plot Functions
#def SlicePlot(field_list, center=(0.5, 0.5, 0.5), width=(2, 'kpc')):
#    plot = yt.SlicePlot(ds, 2, field_list, center, width)
#    plot.save()

#    return 

#def ProjPlot(field_list, center=(0.5, 0.5, 0.5), width=(2, "kpc")):
#    plot = yt.ProjectionPlot(ds, 2, field_list, center, width, weight_field=("gas","density"))
#    plot.save()

#    return

def AddSpeciesField(field_name, atomic_mass):
    def _field_name(field, data): return data[str(field_name + "_density")] / ["Density"] / SolarMetalFractionByMass / atomic_mass
    ds.add_field(("gas", str(field_name)), function=_field_name)

work_dir = "./"
for it_number in range(number_min, number_max+1):
    number  = '%04d' % it_number
    fn = "DD" + number + "/output_" + number


    print(fn)
    ds = yt.load(fn)
    all_data = ds.all_data()
    ds.print_stats()
    print('The current redshift is ' + str(ds.current_redshift))

    print('add field? (Y/N) ')
    add_field = input('')
    if add_field == 'Y' or add_field == 'y' or add_field == 'yes':
        print("what is the species field? ")
        NewSpecies = str(input(''))
        Atomic_Mass = int(input("what is the species atomic mass? "))
        def _NewSpecies(field, data): return data[str(NewSpecies + "_Density")] / data["Density"] / SolarMetalFractionByMass / Atomic_Mass
        ds.add_field(("gas", str(NewSpecies)), function=_NewSpecies, sampling_type="cell", units="1")
    elif add_field == 'N' or add_field == 'n' or len(add_field) == 0:
        pass

    def _cell_size(field, data): return data["cell_volume"]**(1.0/3.0)
    ds.add_field(("gas", "cell_size"), function=_cell_size, sampling_type="cell", units="code_length")
    def _Hydrogen_number_density(field, data): return HydrogenFractionByMass * data["density"] /mp
    ds.add_field(("gas", "Hydrogen_number_density"), function=_Hydrogen_number_density, sampling_type="cell", units="cm**(-3)")
    def _dark_matter_number_density(field, data): return HydrogenFractionByMass * data["dark_matter_density"] /mp
    ds.add_field(("gas", "dark_matter_number_density"), function=_dark_matter_number_density, sampling_type="cell", units="cm**(-3)")
    def _Compressional_heating_rate(field, data): return data["pressure"] * data["velocity_divergence_absolute"]
    ds.add_field(("gas", "Compressional_heating_rate"), function=_Compressional_heating_rate, sampling_type="cell", units="erg/cm**3/s")

    def _Zmet(field, data): return data["SN_Colour"] / data["Density"] / SolarMetalFractionByMass
    ds.add_field(("gas", "Zmet"), function=_Zmet, sampling_type="cell", units="1")
    def _y_H2I(field, data): return data["H2I_Density"] / data["Density"] /HydrogenFractionByMass/2.0
    ds.add_field(("gas", "y_H2I"), function=_y_H2I, sampling_type="cell", units="1")
    def _y_HI(field, data): return data["HI_Density"] / data["Density"] /HydrogenFractionByMass/1.0
    ds.add_field(("gas", "y_HI"), function=_y_HI, sampling_type="cell", units="1")
    def _y_HII(field, data): return data["HII_Density"] / data["Density"] /HydrogenFractionByMass/1.0 
    ds.add_field(("gas", "y_HII"), function=_y_HII, sampling_type="cell", units="1")
    def _y_CII(field, data): return data["CII_Density"] / data["Density"] /HydrogenFractionByMass/12.0 
    ds.add_field(("gas", "y_CII"), function=_y_CII, sampling_type="cell", units="1")
    def _y_OII(field, data): return data["OII_Density"] / data["Density"] /HydrogenFractionByMass/16.0 
    ds.add_field(("gas", "y_OII"), function=_y_OII, sampling_type="cell", units="1")
    def _y_SiI(field, data): return data["SiI_Density"] / data["Density"] /HydrogenFractionByMass/32.0
    ds.add_field(("gas", "y_SiI"), function=_y_SiI, sampling_type="cell", units="1")
    def _y_SiOI(field, data): return data["SiOI_Density"] / data["Density"] /HydrogenFractionByMass/48.0
    ds.add_field(("gas", "y_SiOI"), function=_y_SiOI, sampling_type="cell", units="1")
    def _y_Mg2SiO4(field, data): return data["Mg2SiO4_Density"] / data["Density"] /HydrogenFractionByMass/140.0
    ds.add_field(("gas", "y_Mg2SiO4"), function=_y_Mg2SiO4, sampling_type="cell", units="1")
    def _y_AC(field, data): return data["AC_Density"] / data["Density"] /HydrogenFractionByMass/12.0
    ds.add_field(("gas", "y_AC"), function=_y_AC, sampling_type="cell", units="1")
    def _y_H2O(field, data): return data["H2O_Density"] / data["Density"] /HydrogenFractionByMass/16.0
    ds.add_field(("gas", "y_H2O"), function=_y_H2O, sampling_type="cell", units="1")

    def _ThermalEnergy(field, data): return data["GasEnergy"] * data["cell_mass"]
    ds.add_field(("gas", "ThermalEnergy"), function=_ThermalEnergy, sampling_type="cell", units="erg")
    def _TotEnergy(field, data): return data["TotalEnergy"] * data["cell_mass"]
    ds.add_field(("gas", "TotEnergy"), function=_TotEnergy, sampling_type="cell", units="erg")

    if plot_type == 'Slice':
        print('input plot width ')
        width = float(input('width = '))
        print('input axes unit (mpc, kpc, pc, au...) ')
        axes_unit = str(input('unit = '))
        print('input centre as either as a field or as a coordinate ')
        centre_point = input('input centre point ')
        if centre_point == 'density':
            centre_point = all_data.argmax("density")
        elif centre_point == "Temperature":
            centre_point = all_data.argmax("Temperature")
        elif centre_point == "metallicity":
            centre_point = all_data.argmax("Zmet")
        elif len(centre_point) == 0:
            centre_point = (0.5, 0.5, 0.5) 
        else:
            centre_point = eval(centre_point)

        print(centre_point)

#        centre_point = all_data.argmax("density")
#        print(centre_point)
#        axes_unit = 'kpc'
#        width = (2.0, 'kpc')
#        print("Input field list ")
#        print("Input weight field (e.g. gas) ")
#        field_weight = str(input())
#        print(field_weight)
#        print("Input field (e.g. Temperature) ")
#        field_plot = str(input())
#        print(field_plot)
#        SlicePlot((field_weight, field_plot))
        plot = yt.SlicePlot(ds, 'z', "Zmet", width = width, axes_unit=axes_unit, center=centre_point)
        plot.annotate_grids()
        plot.save("%s/Zmet_%04d.png" % (work_dir, it_number))

        plot = yt.SlicePlot(ds, 'z', "Temperature", width = width, axes_unit=axes_unit, center=centre_point)
        plot.set_log('Temperature', True)
#        plot.set_zlabel(r'$10^{({\rm Temperature})}$')
        plot.save("%s/Temperature_%04d.png" % (work_dir, it_number))

        plot = yt.SlicePlot(ds, 'z', "Hydrogen_number_density", width = width, axes_unit=axes_unit, center=centre_point)
        plot.set_log('Hydrogen_number_density', False)
#        plot.set_zlabel(r'')
        plot.save("%s/H number density_%04d.png" % (work_dir, it_number))

        plot = yt.SlicePlot(ds, 'z', "y_H2O", width = width, axes_unit=axes_unit, center=centre_point)
        plot.annotate_grids()
        plot.save("%s/Water_Abundance_%04d.png" % (work_dir, it_number))

    elif plot_type == 'Proj':
        print('input plot width ')
        width = float(input('width = '))
        print('input axes unit (mpc, kpc, pc, au...) ')
        axes_unit = str(input('unit = '))
        print('input centre as either as a field or as a coordinate ')
        centre_point = input('input centre point ')
        if centre_point == 'density':
            centre_point = all_data.argmax("density")
        elif centre_point == "Temperature":
            centre_point = all_data.argmax("Temperature")
        elif centre_point == "metallicity":
            centre_point = all_data.argmax("Zmet")
        elif len(centre_point) == 0:
            centre_point = (0.5, 0.5, 0.5) 
        else:
            centre_point = eval(centre_point)

        print(centre_point)
        
        print('input weight field (leave blank to default to density) ')
        weight_field = str(input('weight field = '))
        if len(weight_field) <= 0:
            weight_field = "density"
        else:
            weight_field = weight_field

        plot = yt.ProjectionPlot(ds, 'z', "Zmet", width = width, axes_unit=axes_unit, center=centre_point, weight_field=weight_field)
        plot.annotate_grids()
        plot.save("%s/Zmet_%04d.png" % (work_dir, it_number))

        plot = yt.ProjectionPlot(ds, 'z', "Temperature", width = width, axes_unit=axes_unit, center=centre_point, weight_field=weight_field)
        plot.set_log('Temperature', True)
#        plot.set_zlabel(r'$10^{({\rm Temperature})}$')
        plot.save("%s/Temperature_%04d.png" % (work_dir, it_number))

        plot = yt.ProjectionPlot(ds, 'z', "Hydrogen_number_density", width = width, axes_unit=axes_unit, center=centre_point, weight_field=weight_field)
        plot.set_log('Hydrogen_number_density', False)
#        plot.set_zlabel(r'')
        plot.save("%s/H number density_%04d.png" % (work_dir, it_number))

        plot = yt.ProjectionPlot(ds, 'z', "y_H2O", width = width, axes_unit=axes_unit, center=centre_point, weight_field=weight_field)
        plot.annotate_grids()
        plot.save("%s/Water_Abundance_%04d.png" % (work_dir, it_number))

    elif plot_type == 'Phase':
        print('input plot width in kpc ')
        width = int(input('width = '))
        #print('input axes unit (mpc, kpc, pc, au...) ')
        #axes_unit = str(input('unit = '))
        print('input centre as either as a field or as a coordinate ')
        centre_point = input('input centre point ')
        if centre_point == 'density':
            centre_point = all_data.argmax("density")
        elif centre_point == "Temperature":
            centre_point = all_data.argmax("Temperature")
        elif centre_point == "metallicity":
            centre_point = all_data.argmax("Zmet")
        elif len(centre_point) == 0:
            centre_point = (0.5, 0.5, 0.5)
        else:
            centre_point = eval(centre_point)

        print(centre_point)
        
        # Create sphere
        my_sphere = ds.sphere(centre_point, (width, 'kpc')) 

        # Weight field
        print('input weight field (leave blank to default to None) ')
        weight_field = str(input('weight field = '))
        if len(weight_field) <= 0:
            weight_field = None
        else:
            weight_field = weight_field

        # Create the phase plots themselves
        plot = yt.PhasePlot(my_sphere, ("gas", "Hydrogen_number_density"),
        ("Temperature"),
        ("cell_mass"),
        weight_field = weight_field,)
        plot.save("%s/Hydrogen_Number_density_%04d.png" % (work_dir, it_number))



        

        