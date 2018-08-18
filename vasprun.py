#!/usr/bin/env  python
# encoding: utf-8
from lxml import etree
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import numpy as np 

class vasprun:
    """
    parse vasprun.xml and return all useful info to self.values
    Args:
    vasp_file: the path of vasprun.xml
    verbosity: output error msgs or not
    """


    def __init__(self, vasp_file='vasprun.xml', verbosity=0):

        self.error = False
        self.errormsg = ''
        self.values = {}
        try:
            doc = etree.parse(vasp_file)
            doc = doc.getroot()
            self.parse_vaspxml(doc)
            self.get_band_gap()
            N_atom = len(self.values['finalpos']['positions'])
            self.values['calculation']['energy_per_atom'] = \
            self.values['calculation']['energy']/N_atom
        except etree.XMLSyntaxError:
            self.error = True
            self.errormsg = 'corrupted file found'

        if verbosity > 0 and self.error is True:
            print("----------Warning---------------")
            print(self.errormsg)
            print("--------------------------------")
   
    def parse_vaspxml(self, xml_object):
        """
        parse the following tags
        incar
        kpoints
        atominfo - composition, elements, formula
        calculation - eigenvalues, energy, fermi energy, stress, force
        finalpos - lattice, rec_lat, positions
        """

        for child in xml_object.iterchildren():
           self.values.setdefault(child.tag, {})
           if child.tag == "incar":
               self.values[child.tag] = self.parse_i_tag_collection(child)
           elif child.tag == "kpoints":
               self.values[child.tag] = self.parse_kpoints(child)
           elif child.tag == "parameters":
               self.values[child.tag] = self.parse_parameters(child)
           elif child.tag == "atominfo":
               self.values["name_array"] = self.parse_name_array(child)
               self.values["composition"] = self.parse_composition(child)
               self.values["elements"] = self.get_element(self.values["composition"])
               self.values["formula"] = self.get_formula(self.values["composition"])
               self.values["pseudo_potential"], self.values["potcar_symbols"], \
               self.values["valence"], self.values["mass"] = \
                                        self.get_potcar(child)
           elif child.tag == "calculation":
               self.values["calculation"], scf_count = self.parse_calculation(child)
               if self.values['parameters']['response functions']['NELM'] == scf_count:
                   self.error = True
                   self.errormsg = 'SCF is not converged'
           elif child.tag == "structure" and child.attrib.get("name") == "finalpos":
               self.values["finalpos"] = self.parse_finalpos(child)
           elif child.tag not in ("i", "r", "v", "incar", "kpoints", "atominfo", "calculation"):
               self.values[child.tag] = self.parse_vaspxml(child)
           #else:
           #    return 1
           self.dict_clean(self.values)

    @staticmethod
    def dict_clean(dict_del):
        """
        Delete the keys that is {} and None
        Loops recursively over nested dictionaries.
        """
        dict_foo = dict_del.copy()  # Used as iterator to avoid the 'DictionaryHasChanged' error
        for key in dict_foo.keys():
            if isinstance(dict_foo[key], dict):
                vasprun.dict_clean(dict_del[key])

            if dict_foo[key] == {} or dict_foo[key] is None:
                try:
                    del dict_del[key]
                except KeyError:
                    pass

        return dict_del

    def parse_finalpos(self, finalpos):
        """obtain final configuration"""
        d = {}
        for i in finalpos.iter("varray"):
            name = i.attrib.get("name")
            d[name] = self.parse_varray(i)

        return d


    def parse_i_tag_collection(self, itags_collection):
        d = {}
        for info in itags_collection.findall("i"):
            name = info.attrib.get("name")
            type = info.attrib.get("type")
            content = info.text
            d[name] = self.assign_type(type, content)
        return d
   
    @staticmethod
    def parse_varray_pymatgen(elem):
        def _vasprun_float(f):
            """
            Large numbers are often represented as ********* in the vasprun.
            This function parses these values as np.nan
            """
            try:
                return float(f)
            except ValueError as e:
                f = f.strip()
                if f == '*' * len(f):
                    warnings.warn('Float overflow (*******) encountered in vasprun')
                    return np.nan
                raise e
        if elem.get("type", None) == 'logical':
            m = [[True if i=='T' else False for i in v.text.split()] for v in elem]
        else:
            m = [[_vasprun_float(i) for i in v.text.split()] for v in elem]

        return m       

    @staticmethod
    def parse_varray(varray):
        if varray.get("type") == 'int':
            m = [[int(number) for number in v.text.split()] for v in varray.findall("v")]
        else:
            m = [[float(number) for number in v.text.split()] for v in varray.findall("v")]
        return m

    @staticmethod
    def parse_array(array):
        array_dictionary = {}
        values = []
        dimension_list = {}
        field_list = []

        for dimension in array.findall("dimension"):
            dimension_list[dimension.attrib.get("dim")] = dimension.text

        for field in array.findall("field"):
            field_list.append(field.text)

        for r in array.iter("r"):
            values.append([float(number) for number in r.text.split()])

        array_dictionary["value"] = values
        array_dictionary['dimensions'] = dimension_list
        array_dictionary['fileds'] = field_list

        return array_dictionary

    @staticmethod
    def assign_type(type, content):
        if type == "logical":
            if type in ('T', 'True', 'true'):
                return True
            elif type in ('F', 'False', 'false'):
                return False
            else:
                Warning("logical text " + content + " not T, True, true, F, False, false, set to False")
            return False
        elif type == "int":
            return int(content) if len(content.split()) == 1 else [int(number) for number in content.split()]
        elif type == "string":
            return content
        elif type == None:
            return float(content) if len(content.split()) == 1 else [float(number) for number in content.split()]
        else:
            Warning("New type: " + type + ", set to string")
        return content

    @staticmethod
    def parse_composition(atom_info):
        atom_names = {}
        for set in atom_info.findall("array"):

            if set.attrib.get("name") == "atoms":
                for rc in set.iter("rc"):
                    atom_name = rc.find("c").text.replace(" ", '')
                    if atom_name in atom_names:
                        atom_names[atom_name] += 1
                    else:
                        atom_names[atom_name] = 1

                break
        return atom_names

    @staticmethod
    def get_element(atom_names_dictionary):
        elements = []
        for atom_name in atom_names_dictionary:
            elements.append(atom_name.replace(" ", ""))
        return elements

    @staticmethod
    def get_formula(atom_names_dictionary):
        formula = ''
        for atom_name in atom_names_dictionary:
            formula += atom_name.replace(' ', '') + str(atom_names_dictionary[atom_name])
        return formula

    def get_potcar(self, child):
        #{'labels': ['O', 'Sr_sv'], 'pot_type': 'paw', 'functional': 'pbe'}
        #['PAW_PBE', 'N', '08Apr2002']
        pseudo = {'labels': [], 'pot_type':[], 'functional': []}
        potcar_symbol = []
        valence = []
        mass = []
        for i in child.iterchildren():
            if i.tag == "array" and i.attrib.get("name") == 'atomtypes':
               ll = list(i.iter('c'))
               for i in range(3,len(ll),5):
                   valence.append(float(ll[i].text))

               for i in range(2,len(ll),5):
                   mass.append(float(ll[i].text))

               for i in range(4,len(ll),5):
                   text = ll[i].text.split()
                   label = text[1]
                   pot = text[0].split('_')[0]
                   try:
                       xc =  text[0].split('_')[1]
                   except:
                       xc =  'unknown'
                   pseudo['labels'].append(label) 
                   pseudo['pot_type'].append(pot) 
                   pseudo['functional'].append(xc) 
                   potcar_symbol.append(xc + ' ' + label)
        #print(type(valence[0]))          
        return pseudo, potcar_symbol, valence, mass

    @staticmethod
    def parse_name_array(atominfo):
        atom_names = []
        for array in atominfo.findall("array"):
            if array.attrib["name"] == "atoms":
                atom_names = [rc.find("c").text.strip() for rc in array.find("set")]

        if atom_names == []:
            ValueError("No atomname found in file")

        return atom_names

#    def parse_eigenvalue(self, eigenvalue):
#        eigenvalue = eigenvalue.find("array")
#        eigenvalues = self.parse_array(eigenvalue)
#        return eigenvalues
    def parse_parameters(self, child):
        parameters = {}
        for i in child:
            if i.tag == "separator":
                name = i.attrib.get("name")
                d = self.parse_i_tag_collection(i)
                parameters[name]=d
                for ii in i:
                    if ii.tag == "separator":
                        name2 = ii.attrib.get("name")
                        d2 = self.parse_i_tag_collection(ii)
                        parameters[name][name2]=d2
        return parameters

    def parse_eigenvalue(self, eigenvalue):
        eigenvalues = []
        for s in eigenvalue.find("array").find("set").findall("set"):
            for ss in s.findall("set"):
                eigenvalues.append(self.parse_varray_pymatgen(ss))
        return eigenvalues

    def parse_calculation(self, calculation):
        stress = []
        force = []
        efermi = 0.0
        eigenvalues = []
        energy = 0.0
        scf_count = 0
        for i in calculation.iterchildren():
            if i.attrib.get("name") == "stress":
                stress = self.parse_varray(i)
            elif i.attrib.get("name") == "forces":
                force = self.parse_varray(i)
            elif i.tag == "dos":
                for j in i.findall("i"):
                    if j.attrib.get("name") == "efermi":
                        efermi = float(j.text)
                        break
            elif i.tag == "eigenvalues":
                eigenvalues = self.parse_eigenvalue(i)
            elif i.tag == "scstep":
                for j in i.iterchildren():
                    if j.tag == 'energy':
                        for e in j.findall("i"):
                            if e.attrib.get("name") == "e_fr_energy":
                               scf_count +=1

            elif i.tag == "energy":
                for e in i.findall("i"):
                    if e.attrib.get("name") == "e_fr_energy":
                        energy = float(e.text)
                    else:
                        Warning("No e_fr_energy found in <calculation><energy> tag, energy set to 0.0")

        calculation = {}
        calculation["stress"] = stress
        calculation["efermi"] = efermi
        calculation["force"] = force
        calculation["eigenvalues"] = eigenvalues
        calculation["energy"] = energy
        return calculation, scf_count

    def parse_kpoints(self, kpoints):
        kpoints_dict = {'list':[], 'weights':[]}
        for va in kpoints.findall("varray"):
            name = va.attrib["name"]
            if name == "kpointlist":
                kpoints_dict['list'] = self.parse_varray(va)
            elif name == "weights":
                kpoints_dict['weights'] = self.parse_varray(va)
        return kpoints_dict
    
    def get_bands(self):
        """
        Function for computing the valence band index from the count of electrons
        Args: 
            None
        Returns:
            bands: an integer number
            occupy: bool number 
        """
        valence = self.values['valence']
        composition = self.values['composition']
        total = int(self.values['parameters']['electronic']['NELECT'])
        
        if self.values['parameters']['electronic']['electronic spin']['LSORBIT'] == 'T':
           fac = 1
        else:
           fac = 2

        if total%2 == 0:
           IBAND = int(total/fac)
           occupy = True
        else:
           IBAND = int(total/fac) + 1
           occupy = False

        self.values["bands"] = IBAND
        self.values["occupy"] = occupy

    @staticmethod
    def get_cbm(kpoints, efermi, eigens, IBAND):

        ind = np.argmin(eigens[:, IBAND, 0]) 
        pos = kpoints[ind]
        value = eigens[ind, IBAND, 0] - efermi
        return {'kpoint': pos, 'value':value}

    @staticmethod
    def get_vbm(kpoints, efermi, eigens, IBAND):
        
        ind = np.argmax(eigens[:, IBAND, 0]) 
        pos = kpoints[ind]
        value = eigens[ind, IBAND, 0] - efermi

        return {'kpoint': pos, 'value':value}

    def get_band_gap(self):
        self.get_bands()
        IBAND = self.values['bands']
        occupy = self.values['occupy']
        self.values['metal'] = False
        self.values['gap'] = None
        self.values['cbm'] = None
        self.values['vbm'] = None

        if occupy is True:
           efermi = self.values["calculation"]["efermi"] 
           eigens = np.array(self.values['calculation']['eigenvalues'])
           kpoints= np.array(self.values['kpoints']['list'])
           if np.shape(eigens)[0] > np.shape(kpoints)[0]:
              kpoints = np.tile(kpoints, [2,1]) 

           cbm = self.get_cbm(kpoints, efermi, eigens, IBAND)
           vbm = self.get_vbm(kpoints, efermi, eigens, IBAND-1)
           self.values['gap'] = cbm['value'] - vbm['value']
           self.values['cbm'] = cbm
           self.values['vbm'] = vbm
           if self.values['gap'] < 0:
              self.values['metal'] = True
              self.values['gap'] = 0
        else:
           self.values['metal'] = True
           self.values['gap'] = 0

    def eigenvalues_by_band(self, band=0):
        efermi = self.values["calculation"]["efermi"]
        eigens = np.array(self.values['calculation']['eigenvalues'])
        return eigens[:,band,0] - efermi
       
       
    def show_eigenvalues_by_band(self, bands=[0]):
        kpts = self.values['kpoints']['list']
        col_name =  {'K-points': kpts}
        for i, band in enumerate(bands):
            eigen = self.eigenvalues_by_band(band)
            name = 'band'+str(i)
            col_name[name] = eigen
        df = pd.DataFrame(col_name)
        print(tabulate(df, headers='keys', tablefmt='psql'))
 
        
    def export_incar(self, filename=None):
        """export incar"""
        contents = []
        for key in self.values['incar'].keys():
            content = key + ' = ' + str(self.values['incar'][key]) 
            if filename is None:
               print(content)
            else: 
               content += '\n'
               contents.append(str(content))
        if filename is not None:
            with open(filename, 'w') as f:
                 f.writelines(contents)           

    def export_kpoints(self, filename=None):
        """export kpoints"""
        contents = ['KPOINTS\n']
        contents += str(len(self.values['kpoints']['list'])) + '\n'
        contents += ['Cartesian\n']
        for kpt, wt in zip(self.values['kpoints']['list'], self.values['kpoints']['weights']):
            #k1, k2, k3 = kpt
            #print(type(k1), type(wt))
            content = "{:10.4f} {:10.4f} {:10.4f} {:10.4f}".format(kpt[0], kpt[1], kpt[2], wt[0])
            #content = kpt[0] + kpt[1] + kpt[2] + '   ' + wt
            if filename is None:
               print(content)
            else: 
               content += '\n'
               contents.append(str(content))
        if filename is not None:
            with open(filename, 'w') as f:
                 f.writelines(contents)           

    def export_structure(self, filename, fileformat = 'poscar'):
        """export incar"""
        atomNames = self.values["name_array"]
        latt = self.values["finalpos"]["basis"]
        pos = self.values["finalpos"]["positions"]

        struc = Structure(latt, atomNames, pos)
        
        if fileformat == 'poscar':
           Poscar(struc).write_file(filename=filename)
        else:
           CifWriter(struc, symprec=0.01).write_file(filename)

from pprint import pprint
from optparse import OptionParser
import pandas as pd
from tabulate import tabulate


if __name__ == "__main__":
    #------------------------------------------------------------------
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-i", "--incar", dest="incar", metavar='incar file',
                      help="export incar file")
    parser.add_option("-p", "--poscar", dest="poscar",
                      help="export poscar file", metavar="poscar file")
    parser.add_option("-c", "--cif", dest="cif", metavar="cif file",
                      help="export symmetrized cif")
    parser.add_option("-k", "--kpoints", dest="kpoints",
                      help="kpoints file", metavar="kpoints file")
    parser.add_option("-d", "--dosplot", dest="dosplot",
                      help="export dos plot", metavar="dos_plot")
    parser.add_option("-v", "--vasprun", dest="vasprun", default='vasprun.xml',
                      help="path of vasprun.xml file, default: vasprun.xml", metavar="vasprun")
    parser.add_option("-f", "--showforce", dest="force", action='store_true',
                      help="show forces, default: no", metavar="dos_plot")
    parser.add_option("-e", "--eigenvalues", dest="band", action='store_true',
                      help="show eigenvalues in valence/conduction band", metavar="dos_plot")


    (options, args) = parser.parse_args()    
    if options.vasprun is None:
       test = vasprun()
    else:
       test = vasprun(options.vasprun)

  
    # standard output
    if test.values['parameters']['ionic']['NSW'] <=1:
       print('This is a single point calculation')
    #pprint(test.values['kpoints'])

    output = {'formula': None,
              'calculation':['efermi','energy', 'energy_per_atom'],
              'metal': None,
              'gap': None}
    for tag in output.keys():
        if output[tag] is None:
           print(tag, ':  ', test.values[tag])
        else:
           for subtag in output[tag]:
               print(subtag, ':  ', test.values[tag][subtag])
    
    #show VBM and CBM when it is nonmetal
    if test.values['metal'] is False:
       col_name = {'label': ['CBM', 'VBM'],
                   'kpoint':     [test.values['cbm']['kpoint'], test.values['vbm']['kpoint']],
                   'values':     [test.values['cbm']['value'], test.values['vbm']['value']]}
       df = pd.DataFrame(col_name)
       print(tabulate(df, headers='keys', tablefmt='psql'))
   

    col_name = {'valence':    test.values['valence'],
                'labels':     test.values['pseudo_potential']['labels'],
                'functional': test.values['pseudo_potential']['functional']}
    df = pd.DataFrame(col_name)
    print(tabulate(df, headers='keys', tablefmt='psql'))
    
    if options.force: 
        col_name = {'lattice':test.values['finalpos']['basis'],
                   'stress (kbar)': test.values['calculation']['stress']}
        df = pd.DataFrame(col_name)
        #pd.set_option('precision',4)
        print(tabulate(df, headers='keys', tablefmt='psql'))
        col_name = {'atom': test.values['finalpos']['positions'],
                   'force (eV/A)': test.values['calculation']['force']}
        df = pd.DataFrame(col_name)
        print(tabulate(df, headers='keys', tablefmt='psql'))
    
    if options.incar:
        test.export_incar(filename = options.incar)
    elif options.kpoints:
        test.export_kpoints(filename = options.kpoints)
    elif options.poscar:
        test.export_structure(filename = options.poscar)
    elif options.cif:
        test.export_structure(filename = options.cif, fileformat='cif')
    elif options.band:
        vb = test.values['bands']-1
        cb = vb + 1
        test.show_eigenvalues_by_band([vb, cb])
        cbs = test.eigenvalues_by_band(cb)
        vbs = test.eigenvalues_by_band(vb)
        ID = np.argmin(cbs-vbs)
        print("Eigenvalue at CBM: ", min(cbs))
        print("Eigenvalue at VBM: ", max(vbs))
        print("minimum gap at : ", test.values['kpoints']['list'][ID])
        print("CB: ", cbs[ID])
        print("VB: ", vbs[ID])
        print("diff: ", cbs[ID]-vbs[ID])
