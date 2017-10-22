from lxml import etree
from pymatgen.io.cif import CifWriter
from pymatgen import Structure
import numpy as np 

class vasprun:
    """parse vasprun.xml and return all useful info to self.values"""

    def __init__(self, vasp_file='vasprun.xml'):

        self.error = False
        self.values = {}
        try: 
           doc = etree.parse(vasp_file)
           doc = doc.getroot()
           self.parse_vaspxml(doc)
           self.get_band_gap()

        except etree.XMLSyntaxError:
           print('corrupted file found')
           self.error = True
           return None

   
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
               self.values["calculation"] = self.parse_calculation(child)
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
                   xc =  text[0].split('_')[1]
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
        return calculation

    def parse_kpoints(self, kpoints):
        for va in kpoints.findall("varray"):
            name = va.attrib["name"]
            if name == "kpointlist":
                kpoints_matrix = self.parse_varray(va)
        return kpoints_matrix
    
    def get_bands(self):
        valence = self.values['valence']
        composition = self.values['composition']
        total = 0
        for atom, val in zip(composition.keys(), valence):
            total += composition[atom]*val
        
        total = int(total)

        if total%2 == 0:
           IBAND = int(total)/2
           occupy = True
        else:
           IBAND = int(total)/2 + 1
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
        IBAND = int(self.values['bands'])
        occupy = self.values['occupy']
        self.values['metal'] = False
        self.values['gap'] = None
        self.values['cbm'] = None
        self.values['vbm'] = None

        if occupy is True:
           efermi = self.values["calculation"]["efermi"] 
           eigens = np.array(self.values['calculation']['eigenvalues'])
           kpoints= np.array(self.values['kpoints'])
           if np.shape(eigens)[0] > np.shape(kpoints)[0]:
              kpoints = np.tile(kpoints, [2,1]) 

           cbm = self.get_cbm(kpoints, efermi, eigens, IBAND)
           vbm = self.get_vbm(kpoints, efermi, eigens, IBAND-1)
           self.values['gap'] = cbm['value'] - vbm['value']
           self.values['cbm'] = cbm
           self.values['vbm'] = vbm
           if self.values['gap'] < 0:
              self.values['metal'] = True
        else:
           self.values['metal'] = True
        
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

    def export_structure(self, filename, fileformat = 'poscar'):
        """export incar"""
        atomNames = self.values["name_array"]
        latt = self.values["finalpos"]["basis"]
        pos = self.values["finalpos"]["positions"]

        struc = Structure(latt, atomNames, pos)
        
        if fileformat == 'poscar':
           struc.to(fmt='poscar', filename=filename)
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
                      help="export incar")
    parser.add_option("-p", "--poscar", dest="poscar",
                      help="export poscar", metavar="poscar file")
    parser.add_option("-c", "--cif", dest="cif", metavar="cif file",
                      help="export symmetrized cif")
    parser.add_option("-k", "--kpoints", dest="kpoints",
                      help="kpoints list", metavar="kpoints file")
    parser.add_option("-d", "--dosplot", dest="dosplot",
                      help="export dos plot", metavar="dos_plot")
    parser.add_option("-v", "--vasprun", dest="vasprun", default='vasprun.xml',
                      help="path of vasprun.xml file, default: vasprun.xml", metavar="vasprun")
    parser.add_option("-f", "--showforce", dest="force",default='no',
                      help="show forces, default: no", metavar="dos_plot")



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
              'calculation':['efermi','energy'],
              'metal': None,
              'gap': None}
    for tag in output.keys():
        if output[tag] is None:
           print(tag, ':  ', test.values[tag])
        else:
           for subtag in output[tag]:
               print(subtag, ':  ', test.values[tag][subtag])

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
  
    if options.force == 'yes':
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
    elif options.poscar:
       test.export_structure(filename = options.poscar)
    elif options.cif:
       test.export_structure(filename = options.cif)

    #pprint(test.values)
