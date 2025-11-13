import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os
os.chdir(os.path.dirname(__file__))

# Class that stores each pucker state as an object with relevant details to graph
class Dihedral:
    def __init__(self, pucker_state, all_SPE, all_substituents, all_angles, failed):
        self.pucker_state = pucker_state
        self.get_substituents(all_substituents, all_angles, all_SPE)
        self.failed = failed

    def get_substituents(self, all_substituents, all_angles, all_SPE):
        properties_dict = {'Substituent' : all_substituents, 'Angle' : all_angles, 'SPE' : all_SPE}
        self.all_properties = pd.DataFrame(properties_dict)
        self.min_SPE = self.all_properties['SPE'].min()
        self.substituents = []
        max_bound = self.all_properties['Substituent'].iloc[-1] + 1

        for i in range(1, max_bound):
            substituent_df = self.all_properties[self.all_properties['Substituent'] == i]
            substituent_df.sort_values(by='Angle')
            self.substituents.append(substituent_df)

        min_idx = self.all_properties['SPE'].idxmin()
        self.min_substituent = self.all_properties['Substituent'].iloc[min_idx]
        self.min_angle = self.all_properties['Angle'].iloc[min_idx]
        self.initial_substituent = self.all_properties['Substituent'].iloc[0]

    def normalize_SPE(self, min_total_SPE):
        self.normalized_min_SPE = self.min_SPE - min_total_SPE

    def check_energy(self):
        print(self.pucker_state)
        print('------------------------------------')
        finished = False

        for substituent in self.substituents:
            if not(substituent.empty):
                if substituent['SPE'].iloc[0] == self.min_SPE or substituent['SPE'].iloc[len(substituent['SPE'])-1] == self.min_SPE:
                    print('Global Energy Minimum Reached! No More Scans Needed!\n\n')
                    finished = True
                    self.scans_finished = True
                    break

        if not finished:
            print('!!!! More dihedral scans needed !!!!')
            print('Lowest Energy Starting Point for Next Scan:')
            print(f'- Starting Substituent: {self.min_substituent}')
            print(f'- Starting Angle: {self.min_angle}\n\n')
            self.scans_finished = False

    def SPE_dihedral_plot(self):
        colors = ["#4A98D3", '#8E7DBE', "#8EDE4C", '#CD5656', "#FF7A1B"]
        fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(10, 5), sharey=True)
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

        for i, substituent in enumerate(self.substituents):
            axs[i].scatter(substituent['Angle'], substituent['SPE'], c=colors[i])
            axs[i].set_title(f'Substituent {i+1}')

        axs[0].set_xlabel('Dihedral Angle (Degrees)')
        fig.supylabel('Electronic Energy (kcal/mol)')
        fig.suptitle(f'{self.pucker_state} Dihedral Plot')
        plt.tight_layout()
        if self.pucker_state == '1e':
            plt.savefig(f'/Users/sbrooks/Documents/pucker_creator/Visuals/dihedral_plot.svg')
        plt.show()

class Pucker:
    def __init__(self, pucker_state, energies, phi, q, theta, failed, ring_atoms, file_path):
        self.pucker_state = pucker_state
        self.energies = {
            'electronic' : energies[0],
            'electronic_ZPE' : energies[1],
            'electronic_thermal_energy' : energies[2],
            'electronic_thermal_enthalpy' : energies[3],
            'electronic_FE' : energies[4]
        }

        if (pucker_state == '1C4' or pucker_state == '4C1') and not 'unconstrained' in file_path:
            self.phi = 180
        else:
            self.phi = phi
        self.q = q
        self.failed = failed
        self.ring_atoms = ring_atoms
        if ring_atoms == 6:
            self.theta = theta

    def normalize_energies(self, min_energies):
        self.energies['electronic_normalized'] = self.energies['electronic'] - min_energies[0]
        self.energies['electronic_ZPE_normalized'] = self.energies['electronic_ZPE'] - min_energies[1]
        self.energies['electronic_thermal_energy_normalized'] = self.energies['electronic_thermal_energy'] - min_energies[2]
        self.energies['electronic_thermal_enthalpy_normalized'] = self.energies['electronic_thermal_enthalpy'] - min_energies[3]
        self.energies['electronic_FE_normalized'] = self.energies['electronic_FE'] - min_energies[4]

class ReadDihedral:
    def __init__(self, name, file_path):
        self.name = name
        self.file_path = os.path.join(file_path, 'outputs')
        self.get_dihedral_SPE()

    def write_failed(self, state, substituent):
        os.system(f'mkdir {self.file_path}/Failed_Structures/{state}')
        os.system(f'mkdir {self.file_path}/Failed_Structures/{state}/{substituent}')
        os.system(f'cp {self.file_path}/../{state}/base_mol_{substituent}.pdb {self.file_path}/Failed_Structures/{state}')
        os.system(f'cp {self.file_path}/../{state}/{substituent}/{state}.com {self.file_path}/Failed_Structures/{state}/{substituent}')

    def get_dihedral_SPE(self):
        state = ''
        all_SPE = []
        all_angles = []
        all_substituents = []
        substituent = 0
        angle = 0
        SPE = 0
        SPE_string = ''
        state = ''
        failed = False

        os.system(f'mkdir {self.file_path}/Failed_Structures')

        # Holds an object for each pucker state
        self.pucker_states = []

        with open(self.file_path + '/all_energies.txt', 'r') as file:
            for i, line in enumerate(file):
                if 'Pucker' in line:
                    if i != 3 and len(all_SPE) != 0:
                        self.pucker_states.append(Dihedral(state, all_SPE, all_substituents, all_angles, failed))
                        all_SPE = []
                        all_substituents = []
                        all_angles = []
                        failed = False
                    state = line.split()[2]
                elif 'Energy' in line:
                    SPE_string = line.split(': ')[1][:-1]
                    if SPE_string == 'N/A':
                        self.write_failed(state, substituent)
                        all_substituents.pop()
                        all_angles.pop()
                        failed = True
                    else:
                        SPE = float(SPE_string)
                        all_SPE.append(SPE)
                elif 'Substituent' in line:
                    substituent = int(line.split()[1])
                    all_substituents.append(substituent)
                elif 'Angle' in line:
                    angle = int(line.split()[1])
                    if angle > 360:
                        angle -= 360
                    all_angles.append(angle)

            if len(all_SPE) != 0:
                self.pucker_states.append(Dihedral(state, all_SPE, all_substituents, all_angles, failed))
            SPE_all_states = [s.min_SPE for s in self.pucker_states]
            min_total_SPE = min(SPE_all_states)

            for s in self.pucker_states:
                s.normalize_SPE(min_total_SPE)

    def check_energy(self):
        for pucker in self.pucker_states:
            pucker.check_energy()
        self.write_output()

    def write_output(self):
        os.system(f'mkdir {self.file_path}/Next_Scan_Starting_Points')
        os.system(f'mkdir {self.file_path}/Finished_Structures')
        all_output = []
        finished_output = []
        for pucker in self.pucker_states:
            newline = {'Pucker State' : pucker.pucker_state, 
                       'Min SPE' : pucker.min_SPE, 
                       'Min Substituent' : pucker.min_substituent, 
                       'Min Angle' : pucker.min_angle, 
                       'Scan Finished' : pucker.scans_finished}
            all_output.append(newline)
            
            if not(pucker.failed):
                if pucker.scans_finished:
                    os.system(f'cp {self.file_path}/xyz_outputs/{pucker.pucker_state}_{pucker.min_substituent}_{pucker.min_angle}.xyz {self.file_path}/Finished_Structures')
                    finished_output.append(newline)
                else:
                    os.system(f'cp {self.file_path}/xyz_outputs/{pucker.pucker_state}_{pucker.min_substituent}_{pucker.min_angle}.xyz {self.file_path}/Next_Scan_Starting_Points')

        all_output_df = pd.DataFrame(all_output)
        all_output_df.to_excel(f'{self.file_path}/output_energies.xlsx', sheet_name=self.name)
        finished_output_df = pd.DataFrame(finished_output)
        finished_output_df.to_excel(f'{self.file_path}/finished_output_energies.xlsx', sheet_name=self.name)

class ReadPucker:
    def __init__(self, name, file_path):
        self.name = name
        self.file_path = os.path.join(file_path, 'outputs')
        self.get_SPE()

    def get_SPE(self):
        state = ''
        phi = 0.0
        q = 0.0

        electronic = 0
        electronic_ZPE = 0
        electronic_thermal_energy = 0
        electronic_thermal_enthalpy = 0
        electronic_FE = 0

        electronic_string = ''

        state = ''
        failed = False
        theta = -1.0

        os.system(f'mkdir {self.file_path}/Failed_Structures')

        # Holds an object for each pucker state
        self.pucker_states = []
        
        with open(self.file_path + '/all_energies.txt', 'r') as file:
            for i, line in enumerate(file):
                if i == 0:
                    ring_atoms = int(line.split()[2])
                if 'Pucker' in line:
                    if i != 4:
                        energies = [electronic * 627.5095, electronic_ZPE * 627.5095, electronic_thermal_energy * 627.5095,
                                    electronic_thermal_enthalpy * 627.5095, electronic_FE * 627.5095]
                        self.pucker_states.append(Pucker(state, energies, phi, q, theta, failed, ring_atoms, self.file_path))
                        failed = False
                    state = line.split()[2]
                elif 'Electronic Energy' in line:
                    electronic_string = line.split(': ')[1][:-1]
                    if electronic_string == 'N/A':
                        failed = True
                    else:
                        electronic = float(electronic_string)
                elif 'zero-point Energies' in line and not failed:
                    electronic_ZPE = float(line.split(': ')[1][:-1])
                elif 'thermal Energies' in line and not failed:
                    electronic_thermal_energy = float(line.split(': ')[1][:-1])
                elif 'thermal Enthalpies' in line and not failed:
                    electronic_thermal_enthalpy = float(line.split(': ')[1][:-1])
                elif 'Free Energies' in line and not failed:
                    electronic_FE = float(line.split(': ')[1][:-1])
                elif 'Phi' in line:
                    phi = float(line.split()[1])
                elif 'Q' in line:
                    q = float(line.split()[1])
                elif ring_atoms == 6 and 'Theta' in line:
                    theta = float(line.split()[1])

            energies = [electronic * 627.5095, electronic_ZPE * 627.5095, electronic_thermal_energy * 627.5095,
                                    electronic_thermal_enthalpy * 627.5095, electronic_FE * 627.5095]
            self.pucker_states.append(Pucker(state, energies, phi, q, theta, failed, ring_atoms, self.file_path))

            min_electronic = min([pucker.energies['electronic'] for pucker in self.pucker_states])
            min_electronic_ZPE = min([pucker.energies['electronic_ZPE'] for pucker in self.pucker_states])
            min_electronic_thermal_energy = min([pucker.energies['electronic_thermal_energy'] for pucker in self.pucker_states])
            min_electronic_thermal_enthalpy = min([pucker.energies['electronic_thermal_enthalpy'] for pucker in self.pucker_states])
            min_electronic_FE = min([pucker.energies['electronic_FE'] for pucker in self.pucker_states])
            min_energies = [min_electronic, min_electronic_ZPE, min_electronic_thermal_energy,
                            min_electronic_thermal_enthalpy, min_electronic_FE]

            for pucker in self.pucker_states:
                pucker.normalize_energies(min_energies)

    def write_output(self):
        os.system(f'mkdir {self.file_path}/Finished_Structures')
        all_output = []
        finished_output = []
        for pucker in self.pucker_states:
            newline = {'Pucker State' : pucker.pucker_state, 
                        'Electronic (E)' : pucker.energies['electronic'],
                        'E Plus ZPE' : pucker.energies['electronic_ZPE'],
                        'E Plus Thermal Energy' : pucker.energies['electronic_thermal_energy'],
                        'E Plus Thermal Enthalpy' : pucker.energies['electronic_thermal_enthalpy'],
                        'E Plus Free Energy' : pucker.energies['electronic_FE'],
                        'E Norm.' : pucker.energies['electronic_normalized'],
                        'E Plus ZPE Norm.' : pucker.energies['electronic_ZPE_normalized'],
                        'E Plus Thermal Energy Norm.' : pucker.energies['electronic_thermal_energy_normalized'],
                        'E Plus Thermal Enthalpy Norm.' : pucker.energies['electronic_thermal_enthalpy_normalized'],
                        'E Plus Free Energy Norm.' : pucker.energies['electronic_FE_normalized'],
                        'Phi' : pucker.phi,
                        'Q' : pucker.q
            }
            if pucker.ring_atoms == 6:
                newline['Theta'] = pucker.theta

            all_output.append(newline)
            
            if not(pucker.failed):
                os.system(f'cp {self.file_path}/xyz_outputs/{pucker.pucker_state}_opt.xyz {self.file_path}/Finished_Structures')
                finished_output.append(newline)
            else:
                os.system(f'cp {self.file_path}/xyz_outputs/{pucker.pucker_state}_opt.xyz {self.file_path}/Failed_Structures')

        finished_output_df = pd.DataFrame(finished_output)
        finished_output_df.to_excel(f'{self.file_path}/finished_output_energies.xlsx', sheet_name=self.name)