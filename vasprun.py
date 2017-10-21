from lxml import etree

class vasprun:
    """parse vasprun.xml and return all useful info to self.values"""

    def __init__(self, vasp_file='static-vasprun.xml'):

        self.error = False
        self.values = {}
        try: 
           doc = etree.parse(vasp_file)
           doc = doc.getroot()
           self.parse_vaspxml(doc)
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
           elif child.tag == "atominfo":
               self.values["name_array"] = self.parse_name_array(child)
               self.values["composition"] = self.parse_composition(child)
               self.values["elements"] = self.get_element(self.values["composition"])
               self.values["formula"] = self.get_formula(self.values["composition"])
               self.values["pseudo_potential"], self.values["potcar_symbols"] = \
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
        """obtain """
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
        for i in child.iterchildren():
            if i.tag == "array" and i.attrib.get("name") == 'atomtypes':
               ll = list(i.iter('c'))
               for i in range(4,10,5):
                   text = ll[i].text.split()
                   label = text[1]
                   pot = text[0].split('_')[0]
                   xc =  text[0].split('_')[1]
                   pseudo['labels'].append(label) 
                   pseudo['pot_type'].append(pot) 
                   pseudo['functional'].append(xc) 
                   potcar_symbol.append(xc + ' ' + label)
        #print(pp)          
        return pseudo, potcar_symbol 

    @staticmethod
    def parse_name_array(atominfo):
        atom_names = []
        for array in atominfo.findall("array"):
            if array.attrib["name"] == "atoms":
                atom_names = [rc.find("c").text.strip() for rc in array.find("set")]

        if atom_names == []:
            ValueError("No atomname found in file")

        return atom_names

    def parse_eigenvalue(self, eigenvalue):
        eigenvalue = eigenvalue.find("array")
        eigenvalues = self.parse_array(eigenvalue)
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


from pprint import pprint
if __name__ == "__main__":
    
    test = vasprun()
    pprint(test.values)
    #for i in test.values:
    #    print(test.values[i])
