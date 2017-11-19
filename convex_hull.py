from os import path, walk, listdir
from scipy.spatial import ConvexHull
from optparse import OptionParser
from vasprun import vasprun
from pprint import pprint
import pandas as pd
import numpy as np
from tabulate import tabulate
import collections
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def get_subdir(a_dir):
    return sorted([name for name in listdir(a_dir)
            if path.isdir(path.join(a_dir, name))])
#-------------------------------------------------------------------
parser = OptionParser()
parser.add_option("-p",  "--pdir", dest='pdir', default='./',
                  help="by parent directory")
parser.add_option("-f", "--filename", dest="file", default='vasprun.xml',
                  help="vasprun.xml")

(options, args) = parser.parse_args()

total_dir =  get_subdir(options.pdir)
entries = []

for sub_dir in total_dir:
    file_path = sub_dir+'/'+options.file
    if path.exists(file_path):
        vasp_calc = vasprun(file_path)
        if vasp_calc.error is False:
            entries.append(vasp_calc.values)

elements = []
for entry in entries:
    for ele in entry['elements']:
        if ele not in elements:
            elements.append(ele)

print('The system contains the following elements: ', elements)
if len(elements) > 2:
    print('Error, we only support the calculation for binary systems')

else:
    ref = [1e+8]*len(elements)
    for i, ele in enumerate(elements):
        for entry in entries:
            ref0 = entry['calculation']['energy_per_atom']
            #print(entry['elements'])
            if len(entry['elements']) == 1 and entry['elements'][0]==ele and ref[i] > ref0:
                ref[i] = ref0
    ref = np.array(ref)
    ref = np.reshape(ref, [1, len(elements)])
    struc = collections.OrderedDict(
        {'formula':[],'e_dft':[],'e_formation':[],
        'e_above_hull':[],'pts':[]}
        )

    for entry in entries:
        N = np.array([0, 0])
        if entry['composition'].get(elements[0]):
            N[0] = entry['composition'][elements[0]]
        if entry['composition'].get(elements[1]):
            N[1] = entry['composition'][elements[1]]
        x0 = N/sum(N)
        x0 = np.reshape(x0, [len(elements),1])
        x = N[1]/sum(N)
        y0 = entry['calculation']['energy_per_atom'] 
        y1 = y0 - np.dot(ref, x0)
        #print(x, ref, y0, y1)
        struc['formula'].append(entry['formula'])
        #struc['space_group'].append(entry['space_group'])
        struc['e_dft'].append(y0)
        struc['e_formation'].append(y1[0,0])
        struc['pts'].append([x,y1[0,0]])
    #print(struc['pts'])
    pts = np.array(struc['pts'])
    hull = ConvexHull(pts)
    
    for i, pt in enumerate(pts):
        x0, y0 = pt[0], pt[1]
        for simplex in hull.simplices:
            x, y = pts[simplex, 0], pts[simplex, 1]
            #print(x, y, x0, y0)
            if sum(y)<0 and max(y)<=0 and x0>=min(x) and x0<=max(x):
                e_tmp = y0- y[0] - (x0-x[0])*(y[1]-y[0])/(x[1]-x[0]) 
                struc['e_above_hull'].append(e_tmp)
                break
    #pprint(struc)
    df = pd.DataFrame(struc, columns=['formula','e_dft','e_formation','e_above_hull'])
    df = df.sort_values('e_above_hull', ascending=[True])
    print(tabulate(df, headers='keys', tablefmt='psql'))
