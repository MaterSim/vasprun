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
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
plt.style.use("bmh")


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

class units:
    am2kg = 1.6605402e-27
    ev2j = 1.60217733e-19
    plank = 6.626075e-34
    c = 2.99792458e+10
    pi = 3.1415926
    proton_mass = 1836
    ev2hartree = 27.211386245988
    a2bohr = 0.529
    ev2cm = 8065.6

class IR:
    def __init__(self, born_chgs, eigenvalues, eigenvectors, mass, vol):
        self.born_chgs = np.array(born_chgs)
        self.modes = eigenvectors
        self.freqs = eigenvalues
        self.Natom = len(born_chgs)
        IRs = []
        epsilons = []
        for mode, freq in zip(self.modes, self.freqs):
            IRx, IRy, IRz = 0, 0, 0
            mode = np.reshape(mode, [self.Natom, 3])
            for i in range(self.Natom):
                for j in range(3):
                    IRx += mode[i, j] * self.born_chgs[i,j,0]
                    IRy += mode[i, j] * self.born_chgs[i,j,1]
                    IRz += mode[i, j] * self.born_chgs[i,j,2]
            IR = IRx**2 + IRy**2 + IRz**2
            IRs.append(IR)
            if abs(IR) > 1e-2 and abs(freq) > 5e-3:
                # IR active, compute the epsilon
                epsilons.append(compute_epsilon_by_modes(mode, freq, self.born_chgs, vol, mass))
            else:
                epsilons.append(np.zeros([3,3]).flatten())
                #print('IR inactive, skip this mode')
        self.IRs = np.array(IRs)
        self.epsilons = epsilons
    def show(self):
        print("\n   Freq(cm-1)    IR Intensity     E_xx         E_yy         E_zz")
        for ir, freq, eps in zip(self.IRs, self.freqs, self.epsilons):
            print("{:12.3f} {:12.3f} {:12.3f} {:12.3f} {:12.3f}".format(freq*units.ev2cm, ir, eps[0], eps[4], eps[8]))
        eps_sum = np.sum(self.epsilons, axis=0)
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[0], eps_sum[1], eps_sum[2]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[3], eps_sum[4], eps_sum[5]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('Total', eps_sum[6], eps_sum[7], eps_sum[8]))


def compute_epsilon_by_modes(mode, freq, z, V, mass):
    """Compute the epsilon for the given mode"""
    #transform all units to hartree
    freq = freq/units.ev2hartree
    V = V/(units.a2bohr**3)  #bohr
    #mass = mass*units.proton_mass

    # compute the mode effiective charge tensors Z* (3 component)
    zt = np.zeros(3)
    for alpha in range(3):
        for beta in range(3):
            for i in range(len(z)):
                zt[alpha] += z[i, alpha, beta] * mode[i, beta] / np.sqrt(mass[i])
                #zt[alpha] += z[i, alpha, beta] * mode[i, beta]
                #zt[alpha] += z[i, alpha, beta] * mode[i, beta] / mass[i]
    epsilon = np.zeros([3,3])
    # compute the epsilon
    for alpha in range(3):
        for beta in range(3):
            epsilon[alpha, beta] = zt[alpha] * zt[beta] / ((freq)**2)
    factor = 4*units.pi/V/units.proton_mass
    epsilon *= factor
    return epsilon.flatten()


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

            #if len(dict_foo[key])==0:
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


    def parse_dynmat(self, dynmat):
        hessian, eigenvalues, eigenvectors = [], [], []
        for va in dynmat.findall("varray"):
            if va.attrib.get("name") == "hessian":
                hessian = self.parse_varray(va)
            elif va.attrib.get("name") == "eigenvectors":
                eigenvectors = self.parse_varray(va)
        factor = np.sqrt(units.ev2j/1e-20/units.am2kg)
        for v in dynmat.findall("v"):
            for i in v.text.split():
                eigenvalues.append(np.sqrt(abs(float(i)))*factor*units.plank/units.ev2j/2/units.pi)
                #eigenvalues.append(np.sqrt(abs(float(i))))
        return hessian, eigenvalues, eigenvectors

    def parse_born_chg(self, charge):
        """obtain the born charge"""
        chg = []
        for info in charge.findall('set'):
            chg.append(self.parse_varray(info))
        return chg

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
        array_dictionary['fields'] = field_list

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
        if dos.find("partial"):
            for s in dos.find("partial").find("array").findall("set"):
                for i, ss in enumerate(s.findall("set")):
                    p = []
                    for sss in ss.findall("set"):
                        p.append(self.parse_varray_pymatgen(sss))
                    p_dos.append(p)

        return t_dos, p_dos

    def parse_projected(self, proj):
        projected = []
        for s in proj.find("array").find("set").findall("set"):
            for ss in s.findall("set"):
                p = []
                for sss in ss.findall("set"):
                    p.append(self.parse_varray_pymatgen(sss))
                projected.append(p)
        
        return projected #[N_kpts, N_bands]

    def parse_calculation(self, calculation):
        stress = []
        force = []
        efermi = 0.0
        eigenvalues = []
        energy = 0.0
        scf_count = 0
        tdos = []
        pdos = []
        born_chgs = []
        hessian = []
        dyn_eigenvalues = [] 
        dyn_eigenvectors = []
        epsilon_ion = []
        epsilon_ = []
        proj = []
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
            elif i.tag == "projected":
                proj = self.parse_projected(i)
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
            elif i.tag == "array" and i.attrib.get("name") == "born_charges":
                born_chgs = self.parse_born_chg(i)
            elif i.tag == "varray" and i.attrib.get("name") == "epsilon_ion":
                epsilon_ion = self.parse_varray(i)

            elif i.tag == "dynmat":
                hessian, dyn_eigenvalues, dyn_eigenvectors = self.parse_dynmat(i)

        calculation = {}
        calculation["stress"] = stress
        calculation["efermi"] = efermi
        calculation["force"] = force
        calculation["eband_eigenvalues"] = eigenvalues
        calculation["energy"] = energy
        calculation["tdos"] = tdos
        calculation["pdos"] = pdos
        calculation["projected"] = proj
        calculation["born_charges"] = born_chgs
        calculation["hessian"] = hessian
        calculation["normal_modes_eigenvalues"] = dyn_eigenvalues
        calculation["normal_modes_eigenvectors"] = dyn_eigenvectors 
        calculation["epsilon_ion"] = epsilon_ion
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
            eigens = np.array(self.values['calculation']['eband_eigenvalues'])
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
        eigens = np.array(self.values['calculation']['eband_eigenvalues'])
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

    def parse_bandpath(self):
        kpts = self.values['kpoints']['list']
        rec_basis = np.array(self.values['finalpos']['rec_basis'])
       
        def inline(kpt, path):
            if len(path) < 2:
                return True
            else:
                v1 = np.array(kpt) - np.array(path[-1])
                v2 = np.array(path[-1]) - np.array(path[-2])
                v1_norm = np.linalg.norm(v1)
                v2_norm = np.linalg.norm(v2)
                if v1_norm < 1e-3:
                    return False
                else:
                    cos = np.dot(v1, v2)/v1_norm/v2_norm
                    if abs(cos-1) < 1e-2:
                        return True
                    else:
                        #print(kpt, path[-2], path[-1], angle, np.dot(v1, v2)/v1_norm/v2_norm)
                        return False
        paths = []
        path = []
        for i, kpt in enumerate(kpts):
            if inline(kpt, path):
                path.append(kpt)
                if i == len(kpts) - 1:
                    paths.append(path)
            else:
                paths.append(path)
                path = []
                path.append(kpt)

        band_points = []
        pointer = 0
        for i, path in enumerate(paths):
            path = np.array(path)
            dist = np.linalg.norm(np.dot(path[0, :] - path[-1, :], rec_basis))
            x = np.linspace(pointer, pointer+dist, len(path[:, 0]))
            if i == 0:
                band_paths = x
            else:
                band_paths = np.hstack((band_paths, x))
            pointer += dist
            band_points.append(pointer)

        self.values['band_paths'] = band_paths
        self.values['band_points'] = band_points

    def plot_band(self, filename='band.png', styles='normal', ylim=[-20, 3], p_max=1.0):
        efermi = self.values["calculation"]["efermi"]
        eigens = np.array(self.values['calculation']['eband_eigenvalues'])
        paths = self.values['band_paths']
        band_pts = self.values['band_points']
        proj = np.array(self.values["calculation"]["projected"]) #[N_kpts, N_band, Ions, 9]
        cm = plt.cm.get_cmap('RdYlBu')
        nkpt, nband, nocc = np.shape(eigens)
        for i in range(nband):
            band = eigens[:, i, 0] - efermi
            p = np.empty([len(paths)])
            for kpt in range(len(paths)):
                p[kpt] = np.sum(proj[kpt, i, :, :])
            plt.plot(paths, band, c='black', lw=1.0)
            if styles == 'projected':
                p[p>p_max] = p_max
                plt.scatter(paths, band, c=p, vmin=0, vmax=p_max, cmap=cm, s=10)

        for pt in band_pts:
            plt.axvline(x=pt, ls='-', color='k', alpha=0.5)

        if styles == 'projected': 
            cbar = plt.colorbar(orientation='horizontal', 
                                ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
                                )
            cbar.set_label('ratio of projected DOS')#, fontsize=12)
        plt.ylabel("Energy (eV)")
        plt.ylim(ylim)
        plt.xlim([0, paths[-1]])
        plt.xticks([])
        plt.savefig(filename)

    def get_dos(self, rows, style='t'):
        mydos = []
        labels = []
        a_array = self.values["name_array"]
        N_atom = len(a_array)
        tdos = np.array(self.values['calculation']['tdos'])
        pdos = np.array(self.values['calculation']['pdos'])
        a, b, c, d = np.shape(pdos)
        pdos = np.reshape(pdos, [b, a, c, d])
        if style == 't':
            for spin in tdos:
                mydos.append(spin[rows, 1])
                labels.append('total')

        elif style in ['s', 'p', 'd']:
            for spin in pdos:
                spd = spin[0, : , :]
                for i in range(1, N_atom):
                    spd += spin[i, :, :]

                if style == 's':
                    mydos.append(spd[rows, 1])
                elif style == 'p':
                    mydos.append(spd[rows, 2] + spd[rows, 3] + spd[rows, 4])
                else:
                    mydos.append(spd[rows, 5] + spd[rows, 6] + spd[rows, 7] + spd[rows, 8] + spd[rows, 9])
                labels.append(style)

        elif style[0] == 'a':
            if style[1].isdigit():
                ids = style[1].split('-')
                start, end = int(ids[0]), int(ids[1]) 
                ids = range(start, end+1)
            else: 
                ele = style[1:]
                ids = []
                for i in range(N_atom):
                    if ele == a_array[i]:
                        ids.append(i)
            for spin in pdos:
                spd = spin[ids[0], :, :]
                for i in ids[1:]:
                    spd += spin[i, :, :] 
                mydos.append(spd[rows, 1] + spd[rows, 2] + spd[rows, 3] + spd[rows, 4] + \
                             spd[rows, 5] + spd[rows, 6] + spd[rows, 7] + spd[rows, 8] + spd[rows, 9])
                labels.append(style[1:])

        if len(labels) == 2:
            labels[0] += '-up'
            labels[1] += '-down'
            mydos[1] *= -1
        return mydos, labels

    def plot_dos(self, filename='dos.png', smear=None, styles='t', xlim=[-3, 3]):
        """export dos"""
        efermi = self.values['calculation']['efermi']
        tdos = np.array(self.values['calculation']['tdos'][0])
        tdos[:, 0] -= efermi
        e = tdos[:, 0]
        rows = (e > xlim[0]) & (e < xlim[1])
        e = e[rows]
        plt_obj = {}
        for option in styles.split('+'):
            if option == 'spd':
                option = ['s', 'p', 'd']
            else:
                option = [option]
            for style in option:
                mydos, labels = self.get_dos(rows, style)
                for data, label in zip(mydos, labels):
                    plt_obj[label] = data

        fig, ax = plt.subplots()
        lines1 = []
        lines2 = []
        labels1 = []
        labels2 = []
        for label in plt_obj.keys():
            e = np.reshape(e, [len(e), 1])
            data = np.reshape(plt_obj[label], [len(e), 1])
            if smear is not None: 
                data = np.hstack((e, data))
                data = smear_data(data, smear)
                data = data[:, 1]
            if label.find('down') > 0:
                lines2 += ax.plot(e, data)
                labels2.append(label)
            else:
                lines1 += ax.plot(e, data)
                labels1.append(label)
        leg1 = ax.legend(lines1, [label for label in labels1], fancybox=True, loc='upper right')
        leg1.get_frame().set_alpha(0.5)
        if len(lines2) > 0:
            from matplotlib.legend import Legend
            leg2 = Legend(ax, lines2, [label for label in labels2], fancybox=True, loc='lower right')
            ax.add_artist(leg2)
            leg2.get_frame().set_alpha(0.5)

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
    parser.add_option("-b", "--bandplot", dest="bandplot", metavar="band_plot", type=str,
                      help="export band plot, options: normal or projected")
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
    parser.add_option("-n", "--figname", dest="figname", type=str, default='fig.png',
                      help="dos/band figure name, default: fig.png", metavar="figname")
    parser.add_option("-l", "--lim", dest="lim", default='-3,3', 
                      help="dos/band plot lim, default: -3,3", metavar="lim")
    parser.add_option("-m", "--max", dest="max", default=0.5, type=float,
                      help="band plot colorbar, default: 0.5", metavar="max")
    parser.add_option("--dyn", dest="dyn", action='store_true',
                      help="dynamic matrix analysis, default: false", metavar="dyn")


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
        lim = options.lim.split(',')
        lim = [float(i) for i in lim]
        test.plot_dos(styles=options.dosplot, filename=options.figname, xlim=lim, smear=options.smear)
    elif options.bandplot:
        lim = options.lim.split(',')
        lim = [float(i) for i in lim]
        test.parse_bandpath()
        test.plot_band(styles=options.bandplot, filename=options.figname, ylim=lim, p_max=options.max)
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
    elif options.dyn:
        chg = test.values['calculation']['born_charges']
        eig = test.values['calculation']['normal_modes_eigenvalues']
        eigv = test.values['calculation']['normal_modes_eigenvectors']
        vol = np.linalg.det(np.array(test.values['finalpos']['basis']))
        mass = []
        for i, ele in enumerate(test.values["composition"].keys()):
            for m in range(test.values["composition"][ele]):
                mass.append(test.values['mass'][i])
        IR(chg, eig, eigv, mass, vol).show()
        modes = []
        for mode in eigv:
            #mode = np.reshape(mode, [int(len(mode)/3), 3])
            modes.append(np.array(mode))
            #for a in mode:
            #    sum += np.sum(a**2)
            #print(sum)
        print(np.sum(modes[14]*modes[10]))
        eps = np.array(test.values['calculation']['epsilon_ion'])
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('DFPT', eps[0,0], eps[0,1], eps[0,2]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('DFPT', eps[1,0], eps[1,1], eps[1,2]))
        print("{:25s} {:12.3f} {:12.3f} {:12.3f}".format('DFPT', eps[2,0], eps[2,1], eps[2,2]))



