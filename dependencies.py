"""
BioPython dependencies and imports management
"""

from Bio import Seq, SeqUtils, SeqIO, Align, AlignIO, Phylo, Entrez, ExPASy, motifs, Restriction

# Optional imports with graceful fallback
try:
    from Bio import Cluster
except ImportError:
    Cluster = None

try:
    from Bio import KEGG
except ImportError:
    KEGG = None

try:
    from Bio import SearchIO
except ImportError:
    SearchIO = None

try:
    from Bio import SwissProt
except ImportError:
    SwissProt = None

try:
    from Bio.SeqUtils import ProtParam
except ImportError:
    ProtParam = None

try:
    from Bio.Data import CodonTable, IUPACData
except ImportError:
    CodonTable = None
    IUPACData = None

try:
    from Bio import PopGen
except ImportError:
    PopGen = None

try:
    from Bio import Pathway
except ImportError:
    Pathway = None

try:
    from Bio import UniGene
except ImportError:
    UniGene = None

try:
    from Bio import HMM
except ImportError:
    HMM = None

# Bio.SubsMat not available in this installation
MatrixInfo = None

try:
    from Bio import File
except ImportError:
    File = None

from Bio.PDB import PDBParser, PDBIO, Select, Superimposer, NeighborSearch, PPBuilder, is_aa

try:
    from Bio.PDB import DSSP
except ImportError:
    DSSP = None

try:
    from Bio.PDB import HSExposure, ResidueDepth
except ImportError:
    HSExposure = None
    ResidueDepth = None

from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram, BasicChromosome, ColorSpiral

try:
    from Bio.Graphics import Comparative, Distribution
except ImportError:
    Comparative = None
    Distribution = None

from Bio.Graphics.GenomeDiagram import Track

# Standard library and third-party imports
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import pandas as pd
from io import StringIO, BytesIO
import json
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import stats
import logomaker


# Export dictionary of available modules
def get_dependencies():
    """Return dictionary of available BioPython modules"""
    return {
        'Seq': Seq,
        'SeqUtils': SeqUtils,
        'SeqIO': SeqIO,
        'Align': Align,
        'AlignIO': AlignIO,
        'Phylo': Phylo,
        'Entrez': Entrez,
        'ExPASy': ExPASy,
        'motifs': motifs,
        'Restriction': Restriction,
        'Cluster': Cluster,
        'KEGG': KEGG,
        'SearchIO': SearchIO,
        'SwissProt': SwissProt,
        'ProtParam': ProtParam,
        'CodonTable': CodonTable,
        'IUPACData': IUPACData,
        'PopGen': PopGen,
        'Pathway': Pathway,
        'UniGene': UniGene,
        'HMM': HMM,
        'MatrixInfo': MatrixInfo,
        'File': File,
        'PDBParser': PDBParser,
        'PDBIO': PDBIO,
        'Select': Select,
        'Superimposer': Superimposer,
        'NeighborSearch': NeighborSearch,
        'PPBuilder': PPBuilder,
        'is_aa': is_aa,
        'DSSP': DSSP,
        'HSExposure': HSExposure,
        'ResidueDepth': ResidueDepth,
        'SeqRecord': SeqRecord,
        'GenomeDiagram': GenomeDiagram,
        'BasicChromosome': BasicChromosome,
        'ColorSpiral': ColorSpiral,
        'Comparative': Comparative,
        'Distribution': Distribution,
        'Track': Track,
        'plt': plt,
        'sns': sns,
        'nx': nx,
        'np': np,
        'pd': pd,
        'StringIO': StringIO,
        'BytesIO': BytesIO,
        'json': json,
        'KMeans': KMeans,
        'DBSCAN': DBSCAN,
        'AgglomerativeClustering': AgglomerativeClustering,
        'dendrogram': dendrogram,
        'linkage': linkage,
        'stats': stats,
        'logomaker': logomaker
    }
