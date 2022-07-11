# This material was prepared as an account of work sponsored by an agency of the 
# United States Government. Neither the United States Government nor the United 
# States Department of Energy, nor Battelle, nor any of their employees, nor any 
# jurisdiction or organization that has cooperated in the development of these 
# materials, makes any warranty, express or implied, or assumes any legal 
# liability or responsibility for the accuracy, completeness, or usefulness or 
# any information, apparatus, product, software, or process disclosed, or 
# represents that its use would not infringe privately owned rights. Reference 
# herein to any specific commercial product, process, or service by trade name, 
# trademark, manufacturer, or otherwise does not necessarily constitute or imply 
# its endorsement, recommendation, or favoring by the United States Government 
# or any agency thereof, or Battelle Memorial Institute. The views and opinions 
# of authors expressed herein do not necessarily state or reflect those of the 
# United States Government or any agency thereof.
#                    PACIFIC NORTHWEST NATIONAL LABORATORY
#                               operated by
#                                BATTELLE
#                                for the
#                      UNITED STATES DEPARTMENT OF ENERGY
#                       under Contract DE-AC05-76RL01830

import os
import os.path as op
import gzip
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
from ovito.modifiers import SelectTypeModifier, DeleteSelectedModifier, PolyhedralTemplateMatchingModifier, CommonNeighborAnalysisModifier
import numpy as np
import networkx as nx
import logging
import tqdm


class MaterialGraph:
    
    def __init__(self, input, output, 
                 atom_assignment_method='PTM', #'PTM','CNA','aCNA'
                 nondefect_lattice="FCC",
                 cutoff=3.462):
        
        self.ovito_PTM_structure = {'OTHER': PolyhedralTemplateMatchingModifier.Type.OTHER,
                                   'FCC': PolyhedralTemplateMatchingModifier.Type.FCC,
                                   'HCP': PolyhedralTemplateMatchingModifier.Type.HCP,
                                   'BCC': PolyhedralTemplateMatchingModifier.Type.BCC,
                                   'ICO': PolyhedralTemplateMatchingModifier.Type.ICO,
                                   'SC': PolyhedralTemplateMatchingModifier.Type.SC,
                                   'CUBIC_DIAMOND': PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND,
                                   'HEX_DIAMOND': PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND,
                                   'GRAPHENE': PolyhedralTemplateMatchingModifier.Type.GRAPHENE}

        self.ovito_CNA_structure = {'OTHER': CommonNeighborAnalysisModifier.Type.OTHER,
                                   'FCC': CommonNeighborAnalysisModifier.Type.FCC,
                                   'HCP': CommonNeighborAnalysisModifier.Type.HCP,
                                   'BCC': CommonNeighborAnalysisModifier.Type.BCC,
                                   'ICO': CommonNeighborAnalysisModifier.Type.ICO}
        
        self.input = input
        self.output_dir = output
        self.nondefect_lattice = nondefect_lattice
        self.cutoff = cutoff
        self.method = atom_assignment_method
        
        if not op.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        
        self.component_edgelist = op.join(self.output_dir, 'edgelist.txt.gz')
        
        if not op.isfile(self.component_edgelist):
            # assign atom types
            logging.info(f'... assigning atom types via {self.method}')
            self.__assign_atoms(method=self.method, nondefect=self.nondefect_lattice)
            
    
    def __assign_atoms(self, method, nondefect):
        """
        Classify atoms as defective or non-defective
        method (str): method to assign atom types based on lattice
        """
        if method == 'PTM':
            assignment_mod = PolyhedralTemplateMatchingModifier()
            self.remove_type = ovito_PTM_structure[nondefect.upper()]
        elif method == 'CNA':
            #FixedCutoff              
            assignment_mod = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.FixedCutoff,
                                                            cutoff=self.cutoff)
            remove_type = self.ovito_CNA_structure[nondefect.upper()]
            
        elif method == 'aCNA':
            #AdaptiveCutoff
            assignment_mod = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
            remove_type = self.ovito_CNA_structure[nondefect.upper()]
        else:
            logging.warning("Only CNA and PTM are currently implemented.")
            return
        
        pipeline = import_file(self.input)

        
        # compute atom assignments
        pipeline.modifiers.append(assignment_mod)
        
        # select nondefect atoms and remove
        pipeline.modifiers.append(SelectTypeModifier(operate_on = "particles",
                                                     property = "Structure Type",
                                                     types = {remove_type}))
        pipeline.modifiers.append(DeleteSelectedModifier())
        
        data = pipeline.compute()

        
        self.__collect_defect_neighbors(data)
        

    def __collect_defect_neighbors(self, data):            
        # find neighbors of all defect atoms within cutoff radius
        finder = CutoffNeighborFinder(self.cutoff, data)
        
        f = gzip.open(self.component_edgelist, "wb")
        for index in tqdm.tqdm(range(data.particles.count)):
            # iterate over neighbors of the current atom
            for neigh in finder.find(index):
                str_to_write = str(index)+" "+str(neigh.index)+" {'weight': "+str(neigh.distance)+"}\n"
                f.write(bytes(str_to_write, encoding='utf8'))
        f.close()

        
    def subgraph(self):
        self.G = nx.read_edgelist(self.component_edgelist)
        components = sorted(nx.connected_components(self.G), key=len, reverse=True)
        self.C_gbs = components[0]
        self.C_k = components[1:]
        self.C_k_sizes = np.array([len(c) for c in self.C_k])

        
    def monovacancy_check(self):
        # check all components with 12 vertices for isomorphic match with template vacancy
        
        # read in template vacancy
        V = nx.read_edgelist("data/monovacancy.edgelist.txt.gz")
        
        self.confirmed_monovac=[]
        for i in np.where(self.C_k_sizes==12)[0]:
            Ck = self.G.subgraph(self.C_k[i])
            if nx.faster_could_be_isomorphic(V, Ck):
                if nx.is_isomorphic(V, Ck):
                    self.confirmed_monovac+=i