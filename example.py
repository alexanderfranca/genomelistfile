# -*- coding: utf-8 -*-

import pprint

from genomelistfile.genomelistfile import *

gm = GenomeListFile('./tests/fixtures/hsa.list')

gm.generate_map_data() 

pprint.pprint(gm.protein_ecs['hsa:54434'])
