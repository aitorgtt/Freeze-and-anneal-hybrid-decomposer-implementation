'''
* Copyright (c) 2023 TECNALIA <esther.villar@tecnalia.com;eneko.osaba@tecnalia.com>
*
* This file is free software: you may copy, redistribute and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3.
*
* This file is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.

    This is Aitor Gómez-Tejedor´s modification
'''

from pathlib import Path
import numpy as np
import random
import os
import dimod
import networkx as nx
from dwave.system import DWaveSampler


def er_generator(n, p, problem_index=0, save=False):
    if save:
        path = Path(f'instances/ER_{n}_{p}_{problem_index}.txt')
        f = open(path, "w")

    linear = {}
    quadratic = {}

    for i in range(n):
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    edges = 0
    for i in range(n):
        for j in range(i + 1, n):
            var = random.random()
            if var < p:  # RESULTING GRAPH DENSITY
                val = 2 * random.random() - 1
                quadratic[(i, j)] = val
                if save:
                    f.write(f"{i} {j} {val}\n")
                edges += 1
    if save:
        f.close()
        with open(path, 'r') as f:
            with open('newfile.txt', 'w') as f2:
                f2.write(f"{n} {edges}\n")
                f2.write(f.read())
        os.remove(path)
        os.rename('newfile.txt', path)

    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm


def er_native_generator(n, problem_index=0, save=False):
    token = "DEV-e689f2c3ebf1c9ded5c1eb78e31859c7fc31f054"
    sampler = DWaveSampler(token=token)
    qubits = sampler.nodelist
    couplers = sampler.edgelist
    if save:
        path = Path(f'instances/ER_{n}_{0}_{problem_index}.txt')
        f = open(path, "w")

    linear = {}
    quadratic = {}
    used_qubits = []
    k = -100
    for i in qubits:
        if k > 0 and k < n:
            val = 2 * random.random() - 1
            linear[i] = val
            if save:
                f.write(f"{i} {val}\n")
            used_qubits.append(i)
        elif k > n:
            break
        k += 1

    edges = 0
    for v in couplers:
        if v[0] in used_qubits and v[1] in used_qubits:
            val = 2 * random.random() - 1
            quadratic[v] = val
            if save:
                f.write(f"{v[0]} {v[1]} {val}\n")
            edges += 1
    if save:
        f.close()
        with open(path, 'r') as f:
            with open('newfile.txt', 'w') as f2:
                f2.write(f"{n} {edges}\n")
                f2.write(f.read())
        os.remove(path)
        os.rename('newfile.txt', path)

    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm

def er_native_generator_connected(n, problem_index=0, save=False):
    token = "DEV-e689f2c3ebf1c9ded5c1eb78e31859c7fc31f054"
    sampler = DWaveSampler(token=token)
    target_nodelist, target_edgelist, target_adjacency = dimod.child_structure_dfs(sampler)

    P = nx.Graph(target_edgelist) #Pegasus

    radius = 1
    Sold = nx.Graph()
    while True:
        Snew = nx.ego_graph(P, 4293, radius = radius)
        if len(Snew) > n:
            difference = n - len(Sold)
            selected = list(Snew.nodes - Sold.nodes)[:difference] + list(Sold.nodes)
            S = Snew.subgraph(selected)
            break
        elif len(Snew) == n:
            S = Snew
            break
        else:
            radius += 1
            Sold = Snew

    if save:
        edge_number = len(S.edges)
        path = Path(f'instances/ER_{n}_{0}_{problem_index}_connected.txt')
        f = open(path, "w")
        f.write(f"{n} {edge_number}\n")

    linear = {}
    quadratic = {}

    for i in S.nodes:
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    for v in S.edges:
        val = 2 * random.random() - 1
        quadratic[v] = val
        if save:
            f.write(f"{v[0]} {v[1]} {val}\n")
    if save:
        f.close()
    
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm



def OLD_lpb_generator(n, problem_index=0, save=False): #LPB graph = linear probability binomial graph
    ### Deprecated: this is the special case of lnb_generator where p=0.5 always.###
    linear = {}
    quadratic = {}
    if save:
        path = Path(f'instances/lpb_{n}_{0.5}_{problem_index}.txt')
        f = open(path, "w")

    probs = np.linspace(1,0,n*(n-1)//2)

    for i in range(n):
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    edges = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            var = random.random()
            if var < probs[k]:  # RESULTING GRAPH DENSITY
                val = 2 * random.random() - 1
                quadratic[(i, j)] = val
                if save:
                    f.write(f"{i} {j} {val}\n")
                edges += 1
            k += 1
    if save:
        f.close()
        with open(path, 'r') as f:
            with open('newfile.txt', 'w') as f2:
                f2.write(f"{n} {edges}\n")
                f2.write(f.read())
        os.remove(path)
        os.rename('newfile.txt', path)

    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm

def lpb_generator(n, p, problem_index=0, save=False): #LPB graph = linear probability binomial graph

    linear = {}
    quadratic = {}
    if save:
        path = Path(f'instances/lpb_{n}_{p}_{problem_index}.txt')
        f = open(path, "w")

    if p > 0.5:
        probs = np.linspace(1,2*p-1,n*(n-1)//2)
    elif p <= 0.5:
        probs = np.linspace(2*p,0,n*(n-1)//2)

    for i in range(n):
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    edges = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            var = random.random()
            if var < probs[k]:  # RESULTING GRAPH DENSITY
                val = 2 * random.random() - 1
                quadratic[(i, j)] = val
                if save:
                    f.write(f"{i} {j} {val}\n")
                edges += 1
            k += 1
    if save:
        f.close()
        with open(path, 'r') as f:
            with open('newfile.txt', 'w') as f2:
                f2.write(f"{n} {edges}\n")
                f2.write(f.read())
        os.remove(path)
        os.rename('newfile.txt', path)
        
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm


def random_order_lpb_generator(n, p, problem_index=0, save=False):

    linear = {}
    quadratic = {}
    if save:
        path = Path(f'instances/random_order_lpb_{n}_{p}_{problem_index}.txt')
        f = open(path, "w")
        
    if p > 0.5:
        probs = np.linspace(1,2*p-1,n*(n-1)//2)
    elif p <= 0.5:
        probs = np.linspace(2*p,0,n*(n-1)//2)
        
    order = list(range(n))
    random.shuffle(order)
    
    for i in order:
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    edges = 0
    k = 0
    for i in order:
        for j in range(i + 1, n):
            var = random.random()
            if var < probs[k]:  # RESULTING GRAPH DENSITY
                val = 2 * random.random() - 1
                quadratic[(i, j)] = val
                if save:
                    f.write(f"{i} {j} {val}\n")
                edges += 1
            k += 1
    if save:
        f.close()
        with open(path, 'r') as f:
            with open('newfile.txt', 'w') as f2:
                f2.write(f"{n} {edges}\n")
                f2.write(f.read())
        os.remove(path)
        os.rename('newfile.txt', path)

    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm



def ba_generator(n, m, k, p, problem_index=0, save=False):
    """ n : size of graph
        m : size of nucleus
        k : degree of the nodes in the outline
        p : density of the nucleus
    """
    
    token = "DEV-e689f2c3ebf1c9ded5c1eb78e31859c7fc31f054"
    sampler = DWaveSampler(token=token)
    target_nodelist, target_edgelist, target_adjacency = dimod.child_structure_dfs(sampler)

    G = nx.barabasi_albert_graph(n, k, initial_graph = nx.erdos_renyi_graph(m, p))

    if save:
        edge_number = len(G.edges)
        path = Path(f'instances/BA_{n}_{m}_{k}_{p}_{problem_index}.txt')
        f = open(path, "w")
        f.write(f"{n} {edge_number}\n")

    linear = {}
    quadratic = {}

    for i in G.nodes:
        val = 2 * random.random() - 1
        linear[i] = val
        if save:
            f.write(f"{i} {val}\n")

    for v in G.edges:
        val = 2 * random.random() - 1
        quadratic[v] = val
        if save:
            f.write(f"{v[0]} {v[1]} {val}\n")
    if save:
        f.close()
    
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, 0.0, 'BINARY')

    return bqm


