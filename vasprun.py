#!/usr/bin/env  python
# encoding: utf-8
from scipy.ndimage.filters import gaussian_filter1d
from lxml import etree
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import numpy as np
import warnings
from pprint import pprint
from optparse import OptionParser
import pandas as pd
from tabulate import tabulate


def smear_data(data, sigma):
    """
    Apply Gaussian smearing to spectrum y value.

    Args:
        sigma: Std dev for Gaussian smear function
    """
    diff = [data[i + 1, 0] - data[i, 0] for i in range(len(data) - 1)]
    avg_x_per_step = np.sum(diff) / len(diff)
    data[:, 1] = gaussian_filter1d(data[:, 1], sigma / avg_x_per_step)
    return data


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
        calculation - eigenvalues, energy, dos, fermi energy, stress, force
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
            # else:
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
            m = [[True if i == 'T' else False for i in v.text.split()] for v in elem]
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
            content = content.replace(" ", "")
            if content in ('T', 'True', 'true'):
                return True
            elif content in ('F', 'False', 'false'):
                return False
            else:
                Warning("logical text " + content + " not T, True, true, F, False, false, set to False")
            return False
        elif type == "int":
            return int(content) if len(content.split()) == 1 else [int(number) for number in content.split()]
        elif type == "string":
            return content
        elif type is None:
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
        # {'labels': ['O', 'Sr_sv'], 'pot_type': 'paw', 'functional': 'pbe'}
        # ['PAW_PBE', 'N', '08Apr2002']
        pseudo = {'labels': [], 'pot_type': [], 'functional': []}
        potcar_symbol = []
        valence = []
        mass = []
        for i in child.iterchildren():
            if i.tag == "array" and i.attrib.get("name") == 'atomtypes':
                ll = list(i.iter('c'))
                for i in range(3, len(ll), 5):
                    valence.append(float(ll[i].text))

                for i in range(2, len(ll), 5):
                    mass.append(float(ll[i].text))

                for i in range(4, len(ll), 5):
                    text = ll[i].text.split()
                    label = text[1]
                    pot = text[0].split('_')[0]
                    try:
                        xc = text[0].split('_')[1]
                    except:
                        xc = 'unknown'
                    pseudo['labels'].append(label)
                    pseudo['pot_type'].append(pot)
                    pseudo['functional'].append(xc)
                    potcar_symbol.append(xc + ' ' + label)
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
                parameters[name] = d
                for ii in i:
                    if ii.tag == "separator":
                        name2 = ii.attrib.get("name")
                        d2 = self.parse_i_tag_collection(ii)
                        parameters[name][name2] = d2
        return parameters

    def parse_eigenvalue(self, eigenvalue):
        eigenvalues = []
        for s in eigenvalue.find("array").find("set").findall("set"):
            for ss in s.findall("set"):
                eigenvalues.append(self.parse_varray_pymatgen(ss))
        return eigenvalues

    def parse_dos(self, dos):
        t_dos = []
        p_dos = []
        for s in dos.find("total").find("array").findall("set"):
            for ss in s.findall("set"):
                t_dos.append(self.parse_varray_pymatgen(ss))
        for s in dos.find("partial").find("array").findall("set"):
            for ss in s.findall("set"):
                for sss in ss.findall("set"):
                    p_dos.append(self.parse_varray_pymatgen(sss))
        return t_dos, p_dos

    def parse_calculation(self, calculation):
        stress = []
        force = []
        efermi = 0.0
        eigenvalues = []
        energy = 0.0
        scf_count = 0
        tdos = []
        pdos = []

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
                tdos, pdos = self.parse_dos(i)
            elif i.tag == "eigenvalues":
                eigenvalues = self.parse_eigenvalue(i)
            elif i.tag == "scstep":
                for j in i.iterchildren():
                    if j.tag == 'energy':
                        for e in j.findall("i"):
                            if e.attrib.get("name") == "e_fr_energy":
                                scf_count += 1

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
        calculation["tdos"] = tdos
        calculation["pdos"] = pdos

        return calculation, scf_count

    def parse_kpoints(self, kpoints):
        kpoints_dict = {'list': [], 'weights': []}
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

        if self.values['parameters']['electronic']['electronic spin']['LSORBIT']:
            fac = 1
        else:
            fac = 2

        if total % 2 == 0:
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
        return {'kpoint': pos, 'value': value}

    @staticmethod
    def get_vbm(kpoints, efermi, eigens, IBAND):
        ind = np.argmax(eigens[:, IBAND, 0])
        pos = kpoints[ind]
        value = eigens[ind, IBAND, 0] - efermi
        return {'kpoint': pos, 'value': value}

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
            kpoints = np.array(self.values['kpoints']['list'])
            if np.shape(eigens)[0] > np.shape(kpoints)[0]:
                kpoints = np.tile(kpoints, [2, 1])

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
        return eigens[:, band, 0] - efermi

    def show_eigenvalues_by_band(self, bands=[0], spin=True):
        kpts = self.values['kpoints']['list']
        col_name = {'K-points': kpts}
        for i, band in enumerate(bands):
            eigen = self.eigenvalues_by_band(band)
            if spin:
                eigens = np.reshape(eigen, [int(len(eigen)/2), 2])
                name1 = 'band' + str(i) + 'up'
                name2 = 'band' + str(i) + 'down'
                col_name[name1] = eigens[:, 0]
                col_name[name2] = eigens[:, 1]
            else:
                name = 'band'+str(i)
                col_name[name] = eigen
        df = pd.DataFrame(col_name)
        print(df)

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
            content = "{:10.4f} {:10.4f} {:10.4f} {:10.4f}".format(kpt[0], kpt[1], kpt[2], wt[0])
            if filename is None:
                print(content)
            else:
                content += '\n'
                contents.append(str(content))
        if filename is not None:
            with open(filename, 'w') as f:
                f.writelines(contents)

    def export_structure(self, filename, fileformat='poscar'):
        """export incar"""
        atomNames = self.values["name_array"]
        latt = self.values["finalpos"]["basis"]
        pos = self.values["finalpos"]["positions"]

        struc = Structure(latt, atomNames, pos)

        if fileformat == 'poscar':
            Poscar(struc).write_file(filename=filename)
        else:
            CifWriter(struc, symprec=0.01).write_file(filename)

    def plot_dos(self, filename='dos.png', smear=None, styles='t', xlim=[-3, 3]):
        """export dos"""
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
        plt.style.use("bmh")

        efermi = self.values['calculation']['efermi']
        tdos = np.array(self.values['calculation']['tdos'][0])
        tdos[:, 0] -= efermi
        e = tdos[:, 0]
        rows = (e > xlim[0]) & (e < xlim[1])
        e = e[rows]
        plt_obj = {}
        for option in styles.split('+'):
            if option == 't':
                if len(self.values['calculation']['tdos']) == 1:
                    tdos = tdos[rows, :]
                    plt_obj['total'] = tdos[:, 1]
                else:
                    tdos1 = np.array(self.values['calculation']['tdos'][0])
                    tdos2 = np.array(self.values['calculation']['tdos'][1])
                    plt_obj['total-up'] = tdos1[rows, 1]
                    plt_obj['total-down'] = -1*tdos2[rows, 1]
            elif option == 'spd':
                pdos = self.values['calculation']['pdos']
                N_atom = len(self.values["name_array"])
                if len(pdos) == N_atom:
                    pdos = np.array(pdos)
                    spd = pdos[0, :, :]
                    for i in range(1, N_atom):
                        spd += pdos[i, :, :]
                    s = spd[rows, 1]
                    p = spd[rows, 2] + spd[rows, 3] + spd[rows, 4]
                    d = spd[rows, 5] + spd[rows, 6] + spd[rows, 7] + spd[rows, 8] + spd[rows, 9]
                    plt_obj['s'] = s
                    plt_obj['p'] = p
                    plt_obj['d'] = d
                else:
                    pdos = np.array(pdos)
                    spd1 = pdos[0, :, :]
                    spd2 = pdos[1, :, :]
                    for i in range(2, 2*N_atom):
                        if i % 2 == 1:
                            spd1 += np.array(pdos[i, :, :])
                        else:
                            spd2 += np.array(pdos[i, :, :])
                    s1 = spd1[rows, 1]
                    p1 = spd1[rows, 2] + spd1[rows, 3] + spd1[rows, 4]
                    d1 = spd1[rows, 5] + spd1[rows, 6] + spd1[rows, 7] + spd1[rows, 8] + spd1[rows, 9]
                    s2 = spd2[rows, 1]
                    p2 = spd2[rows, 2] + spd2[rows, 3] + spd2[rows, 4]
                    d2 = spd2[rows, 5] + spd2[rows, 6] + spd2[rows, 7] + spd2[rows, 8] + spd2[rows, 9]
                    plt_obj['s-up'] = s1
                    plt_obj['p-up'] = p1
                    plt_obj['d-up'] = d1
                    plt_obj['s-down'] = -1*s2
                    plt_obj['p-down'] = -1*p2
                    plt_obj['d-down'] = -1*d2

        fig, ax = plt.subplots()
        lines1 = []
        lines2 = []
        labels1 = []
        labels2 = []
        #ax.axis('equal')
        for label in plt_obj.keys():
            # print(len(e), len(plt_obj[label]), label)
            e = np.reshape(e, [len(e), 1])
            data = np.reshape(plt_obj[label], [len(e), 1])
            if smear is not None: 
                data = np.hstack((e, data))
                data = smear_data(data, smear)
                data = data[:, 1]
            #print(label, label.find('down'))
            if label.find('down') > 0:
                lines2 += ax.plot(e, data)#, label=label)
                labels2.append(label)
            else:
                lines1 += ax.plot(e, data)#, label=label)
                labels1.append(label)
        ax.legend(lines1, [label for label in labels1], loc='upper right')
        if len(lines2) > 0:
            from matplotlib.legend import Legend
            leg = Legend(ax, lines2, [label for label in labels2], loc='lower right')

        ax.add_artist(leg)
        plt.xlabel("Energy (eV)")
        plt.ylabel("DOS")
        plt.xlim(xlim)
        plt.savefig(filename)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--incar", dest="incar", action='store_true',
                      help="export incar file", metavar='incar file')
    parser.add_option("-p", "--poscar", dest="poscar",
                      help="export poscar file", metavar="poscar file")
    parser.add_option("-c", "--cif", dest="cif", metavar="cif file",
                      help="export symmetrized cif")
    parser.add_option("-k", "--kpoints", dest="kpoints", action='store_true',
                      help="kpoints file", metavar="kpoints file")
    parser.add_option("-d", "--dosplot", dest="dosplot", metavar="dos_plot", type=str,
                      help="export dos plot, options: t, spd, a, a-Si, a-1")
    parser.add_option("-v", "--vasprun", dest="vasprun", default='vasprun.xml',
                      help="path of vasprun.xml file, default: vasprun.xml", metavar="vasprun")
    parser.add_option("-f", "--showforce", dest="force", action='store_true',
                      help="show forces, default: no", metavar="dos_plot")
    parser.add_option("-a", "--allparameters", dest="parameters", action='store_true',
                      help="show all parameters", metavar="parameter")
    parser.add_option("-e", "--eigenvalues", dest="band", action='store_true',
                      help="show eigenvalues in valence/conduction band", metavar="dos_plot")
    parser.add_option("-s", "--smear", dest="smear", type='float',
                      help="smearing parameter for dos plot, e.g., 0.1 A", metavar="smearing")
    parser.add_option("-n", "--dosfig", dest="dosfig", type='float',
                      help="dos figure name, e.g., dos.png", metavar="dosfig")

    (options, args) = parser.parse_args()
    if options.vasprun is None:
        test = vasprun()
    else:
        test = vasprun(options.vasprun)

    # standard output
    if test.values['parameters']['ionic']['NSW'] <= 1:
        print('This is a single point calculation')
    # pprint(test.values['kpoints'])

    output = {'formula': None,
              'calculation': ['efermi', 'energy', 'energy_per_atom'],
              'metal': None,
              'gap': None}
    for tag in output.keys():
        if output[tag] is None:
            print(tag, ':  ', test.values[tag])
        else:
            for subtag in output[tag]:
                print(subtag, ':  ', test.values[tag][subtag])

    # show VBM and CBM when it is nonmetal
    if test.values['metal'] is False:
        col_name = {'label': ['CBM', 'VBM'],
                    'kpoint': [test.values['cbm']['kpoint'], test.values['vbm']['kpoint']],
                    'values': [test.values['cbm']['value'], test.values['vbm']['value']]}
        df = pd.DataFrame(col_name)
        print(tabulate(df, headers='keys', tablefmt='psql'))

    col_name = {'valence': test.values['valence'],
                'labels': test.values['pseudo_potential']['labels'],
                'functional': test.values['pseudo_potential']['functional']}
    df = pd.DataFrame(col_name)
    print(tabulate(df, headers='keys', tablefmt='psql'))

    if options.force:
        col_name = {'lattice': test.values['finalpos']['basis'],
                    'stress (kbar)': test.values['calculation']['stress']}
        df = pd.DataFrame(col_name)
        print(tabulate(df, headers='keys', tablefmt='psql'))
        col_name = {'atom': test.values['finalpos']['positions'],
                    'force (eV/A)': test.values['calculation']['force']}
        df = pd.DataFrame(col_name)
        print(tabulate(df, headers='keys', tablefmt='psql'))

    if options.incar:
        test.export_incar(filename=options.incar)
    elif options.kpoints:
        test.export_kpoints(filename=options.kpoints)
    elif options.poscar:
        test.export_structure(filename=options.poscar)
    elif options.cif:
        test.export_structure(filename=options.cif, fileformat='cif')
    elif options.parameters:
        pprint(test.values['parameters'])
    elif options.dosplot:
        test.plot_dos(styles=options.dosplot, smear=options.smear)
    elif options.band:
        vb = test.values['bands']-1
        cb = vb + 1
        test.show_eigenvalues_by_band([vb, cb])
        cbs = test.eigenvalues_by_band(cb)
        vbs = test.eigenvalues_by_band(vb)
        ID = np.argmin(cbs-vbs)
        if len(cbs) == len(test.values['kpoints']['list']):
            print("Eigenvalue at CBM: ", min(cbs))
            print("Eigenvalue at VBM: ", max(vbs))
            print("minimum gap at : ", test.values['kpoints']['list'][ID])
            print("CB: ", cbs[ID])
            print("VB: ", vbs[ID])
            print("diff: ", cbs[ID]-vbs[ID])
        else:
            print("This is spin calculation")
