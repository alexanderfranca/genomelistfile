import re
import sys
import pprint
import glob
import os
import datetime


class GenomeListFile:
    """
    This class deals with the relation between organisms and metabolic pathway maps from KEGG.

    That relation is stored in files like hsa.list, lma.list etc, inside the directory pathway/organisms/

    Inside that directory there's all organisms gziped files, concerning metabolic pathway maps.

    But one file is special: the file $organism.list (where $organism is the three letter organism code from KEGG).

    Inside the $organism.list file we find lines like below (in the example the organism is zmc and the file is zmc.list):

        path:zmc00010   zmc:A265_00092  zmc:gpm_a ko:K01834 ec:5.4.2.11
        path:zmc00010   zmc:A265_00094  zmc:adh_a ko:K13953 ec:1.1.1.1
        path:zmc00010   zmc:A265_00122  zmc:pgi ko:K01810 ec:5.3.1.9
        path:zmc00010   zmc:A265_00420  ko:K01785 ec:5.1.3.3
        path:zmc00010   zmc:A265_00761  zmc:lpd ko:K00382 ec:1.8.1.4
        path:zmc00010   zmc:A265_00763  zmc:pdh_c ko:K00627 ec:2.3.1.12 ec:3.1.99.- ec:2.7.-.- ec:1.1.1.1

    """

    def __init__(self, file_to_parse):

        self.file_to_parse = file_to_parse
        self.protein_ecs   = {}
        self.protein_maps  = {}
        self.map_ecs       = {}
        self.organism_maps = [] 
        self.organism_ec_numbers = []
        self.pathway_relations_directory = None

        self.generate_map_data()

    def ec_numbers_from_entry( self, string=None ):
        """
        Returns the EC numbers found in entry.
        
        That kind of entry is something like: 'path:hsa00010   hsa:2023    hsa:ENO1 ko:K01689 ec:4.2.1.11'

        An entry line can have several EC numbers separated by blank space.

        Args:

            string(str): An entry line.

        Returns:

            (list): List of found EC numbers in the line.
        """
        
        # Some comments, just in case you have to test any weird behavior.
        # string = 'path:hsa00052   hsa:2683    hsa:B4GALT1 ko:K07966 ec:2.4.1.- ec:2.4.1.38 ec:2.4.1.90 ec:2.4.1.22'
        # expected_ecs = [ '2.4.1.-', '2.4.1.38', '2.4.1.90', '2.4.1.22' ]

        find_ec = re.compile('^path:[a-z].*?\sec:(.*)')

        result = find_ec.search( string )

        if result:
            ecs_found = result.group(1)
            ecs_found = ecs_found.replace('ec:', '')
            ecs_found = ecs_found.split(' ')
        else:
            ecs_found = []

        return ecs_found 


    def protein_identification_from_entry( self, string=None ):
        """
        Returns the protein identification found in an entry.
        
        That kind of entry is something like: 'path:hsa00010   hsa:2023    hsa:ENO1 ko:K01689 ec:4.2.1.11'

        An entry line can have only a single protein identification (hsa:2023 for example).

        Args:

            string(str): An entry line.

        Returns:

            (list): The protein identification.
        """
        
        # Some comments, just in case you have to test any weird behavior.
        # string = 'path:hsa00052   hsa:2683    hsa:B4GALT1 ko:K07966 ec:2.4.1.- ec:2.4.1.38 ec:2.4.1.90 ec:2.4.1.22'
        # expected_ecs = [ '2.4.1.-', '2.4.1.38', '2.4.1.90', '2.4.1.22' ]
        # expected_protein_identification = 'hsa:2683' 

        find_protein_id = re.compile('path:.*?\s([^\s].*?)\s')

        result = find_protein_id.search( string )

        if result:
            protein_identification = result.group(1)
            protein_identification = protein_identification.replace(' ', '')
        else:
            protein_identification = None 

        return protein_identification 


    def metabolic_map_number( self, string=None ):
        """
        Returns the metabolic map number from a string.

        That kind of string is something like: 'path:hsa00010   hsa:2023    hsa:ENO1 ko:K01689 ec:4.2.1.11'
        """

        re_map  = re.compile("^path:[a-z]{1,}([0-9]{1,})\s.*$")
        result = re_map.search( string )

        return result.group(1)


    def generate_map_data( self ):
        """
        Read $organism.list and fills some important dictionaries: protein_ecs, protein_maps, map_ecs and organism_maps.

        Args:
            list_file(str): Full path for the $org.list file.
        Returns:

            (void): Only create and fills class dictionaries. 
        """

        # ----------------------------------------------------------------------------------------------------------------------
        #
        # Not so important... but really important: protein/genes identification will be converted to lower case. Why?
        # We never know the case of the annotatoins.
        # 
        # ----------------------------------------------------------------------------------------------------------------------
        self.organism_maps = []
        self.protein_ecs   = {}
        self.protein_maps  = {}
        self.map_ecs       = {}
        self.organism_maps = [] 
        self.all_ec_numbers = []

        with open(self.file_to_parse, 'r') as f:

            for line in f:
                line = line.rstrip('\r\n')

                map_number = self.metabolic_map_number( line )
                protein_identification = self.protein_identification_from_entry( line )
                protein_identification = protein_identification.lower()

                if not protein_identification in self.protein_ecs: 
                    self.protein_ecs[ protein_identification ] = [] 

                if not protein_identification in self.protein_maps: 
                    self.protein_maps[ protein_identification ] = [] 

                if not map_number in self.map_ecs:
                    self.map_ecs[ map_number ] = []

                ec_numbers = self.ec_numbers_from_entry( line )

                # Will be used as 'all ec numbers from the organism'
                self.all_ec_numbers.extend( ec_numbers )

                self.map_ecs[ map_number ].extend( ec_numbers )
                self.map_ecs[ map_number ] = set( self.map_ecs[ map_number ] )
                self.map_ecs[ map_number ] = list( self.map_ecs[ map_number ] )

                self.organism_maps.append( map_number )
                
                self.protein_maps[ protein_identification ].append( map_number )
                self.protein_maps[ protein_identification ] = set( self.protein_maps[ protein_identification ] )
                self.protein_maps[ protein_identification ] = list( self.protein_maps[ protein_identification ] )

                self.protein_ecs[ protein_identification ].extend( ec_numbers ) 
                self.protein_ecs[ protein_identification ] = set( self.protein_ecs[ protein_identification ] )
                self.protein_ecs[ protein_identification ] = list( self.protein_ecs[ protein_identification ] )

        self.organism_maps = set( self.organism_maps )
        self.organism_maps = list( self.organism_maps )

        self.all_ec_numbers = set(self.all_ec_numbers)
        self.all_ec_numbers = list(self.all_ec_numbers)

        self.organism_ec_numbers = self.all_ec_numbers

        return self.all_ec_numbers


    def _protein_ecs( self ):
        """
        Returns a dictionary containing proteins and its EC numbers.

        Returns:

            (dict): Dictionary with all protein EC numbers.
        """

        return self.protein_ecs


    def ec_by_protein_identification( self, gene=None ):
        """
        Returns EC numbers by protein identification.

        Returns:

            (list): EC numbers.
        """

        if not gene in self.protein_ecs:
            ecs = [] 
        else:
            ecs = self.protein_ecs[ gene ]

        return ecs


    def _protein_maps( self ):
        """
        Returns a dictionary containing proteins and its metabolic pathway maps.

        Returns:

            (dict): Dictionary with all metabolic maps number.
        """

        return self.protein_maps


    def _map_ecs( self ):
        """
        Returns a dictionary containing metabolic pathway maps number and its EC numbers.

        Returns:

            (dict): Dictionary with all metabolic maps number and EC numbers.
        """

        return self.map_ecs


    def _organism_maps( self ):
        """
        Returns a list containing metabolic pathway maps number.

        Returns:

            (dict): Llist with all metabolic maps number.
        """

        return self.organism_maps


    def all_ec_numbers( self ):
        """
        Returns all EC numbers found in the file.

        Since the file is organisms related, that's ok to think that this method returns all EC numbers by organism.

        Returns:

            (list): All EC numbers in the file.
        """

        return self.organism_ec_numbers 


    def all_incomplete_ec_numbers( self ):
        """
        Returns all incomplete EC numbers found in the file.

        Returns:

            (list): Incomplete EC numbers
        """

        incomplete_ecs = []

        re_incomplete_ec = re.compile('-')

        ecs = self.all_ec_numbers

        for ec in ecs:
            if re_incomplete_ec.search( ec ):
                incomplete_ecs.append( ec )

        return incomplete_ecs


    def all_complete_ec_numbers( self ):
        """
        Returns all incomplete EC numbers found in the file.

        Returns:

            (list): Incomplete EC numbers
        """

        complete_ecs = []

        re_incomplete_ec = re.compile('-')

        ecs = self.all_ec_numbers()

        for ec in ecs:
            if not re_incomplete_ec.search( ec ):
                complete_ecs.append( ec )

        return complete_ecs


    # TODO: test, comment
    def _organism_ec_numbers( self ):

        return self.organism_ec_numbers


##    def all_map_data(self):
##
###        organism_maps = []
###        protein_ecs   = {}
###        protein_maps  = {}
###        map_ecs       = {}
###        organism_maps = [] 
###        all_ec_numbers = []
###
##
##        organisms_and_ecs = {}
##        organisms_and_maps = {}
##        all_ec_numbers = []
##        ecs_and_organisms = {}
##
##        for list_file in glob.glob(self.list_files_directory + '/*.list'):
##
##            genome_code = os.path.basename(list_file )
##            genome_code = re.sub('.list','', genome_code)
##
##            self.generate_map_data(list_file)
##
##            organisms_and_ecs[genome_code] = self._organism_ec_numbers()
##            organisms_and_maps[genome_code] = self._organism_maps()
##            all_ec_numbers.extend(self.all_ec_numbers())
##
##        for organism, ecs in organisms_and_ecs.iteritems():
##            
##            for ec in ecs:
##                if not ec in ecs_and_organisms:
##                    ecs_and_organisms[ec] = []
##
##                ecs_and_organisms[ec].append(organism)
##
##        ecs_and_organisms = {}
##
##        all_ec_numbers = set(all_ec_numbers)
##        all_ec_numbers = list(all_ec_numbers)
##
##        return { 'all_ec_numbers': all_ec_numbers, 'ecs_and_genomes': ecs_and_organisms, 'genome_and_ecs': organisms_and_ecs, 'genome_and_maps': organisms_and_maps }

