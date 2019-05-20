import sys
import os
import unittest
from genomelistfile.genomelistfile import *
import re


class TestGenomeListFile( unittest.TestCase ):

    def setUp( self ):
        self.gm = GenomeListFile(file_to_parse='./tests/fixtures/hsa.list')

    def test_metabolic_map_number( self ):

        result =  self.gm.metabolic_map_number('path:hsa00010   hsa:2023    hsa:ENO1 ko:K01689 ec:4.2.1.11')
        self.assertEqual( result, '00010' )

    def test_generate_map_data( self ):

        self.assertTrue( type( self.gm.protein_ecs ) is dict )
        self.assertTrue( type( self.gm.protein_maps ) is dict )
        self.assertTrue( type( self.gm.map_ecs ) is dict )
        self.assertTrue( type( self.gm.organism_maps) is list )

        expected = [ '3.2.1.183', '2.7.1.60']
        result = self.gm.protein_ecs['hsa:10020'] 

        self.assertTrue( set(expected) == set(result) ) 

        expected = ['00232', '00983', '05204', '01100']
        result = self.gm.protein_maps['hsa:10']

        self.assertTrue( set(expected) == set(result) ) 

        expected = ['2.4.1.17', '1.2.1.31', '1.1.1.22', '1.2.1.47', '1.13.99.1', '1.2.1.8', '1.2.1.3', '3.1.1.17']
        result = self.gm.map_ecs['00053']

        self.assertTrue( set(expected) == set(result) )

        self.assertEqual( len(self.gm.organism_maps), 318 )

    def test_maps_data( self ):

        self.assertTrue( type( self.gm.organism_maps) is list )


    def test_ec_numbers_from_entry( self ):

        #string = 'path:hsa00010   hsa:501 hsa:ALDH7A1 ko:K14085 ec:1.2.1.3 ec:1.2.1.8 ec:1.2.1.31'

        string = 'path:hsa00052   hsa:2683    hsa:B4GALT1 ko:K07966 ec:2.4.1.- ec:2.4.1.38 ec:2.4.1.90 ec:2.4.1.22'

        expected_ecs = [ '2.4.1.-', '2.4.1.38', '2.4.1.90', '2.4.1.22' ]

        ecs = self.gm.ec_numbers_from_entry( string )

        self.assertTrue( type( ecs ) is list )
        self.assertEqual( ecs, expected_ecs )


    def test_protein_identification_from_entry( self ):

        string = 'path:hsa00052   hsa:2683    hsa:B4GALT1 ko:K07966 ec:2.4.1.- ec:2.4.1.38 ec:2.4.1.90 ec:2.4.1.22'

        expected_protein_id = 'hsa:2683' 

        protein_id = self.gm.protein_identification_from_entry( string )

        self.assertTrue( type( protein_id ) is str )
        self.assertEqual( protein_id, expected_protein_id )


    def test_protein_ecs( self ):

        data = self.gm._protein_ecs()

        self.assertTrue( type( data ) is dict )
        self.assertTrue( len( data ) > 1 )

        expected = [ '3.2.1.183', '2.7.1.60']
        result = self.gm.protein_ecs['hsa:10020'] 

        self.assertTrue( set(expected) == set(result) ) 

    def test_protein_maps( self ):

        data = self.gm._protein_maps()

        self.assertTrue( type( data ) is dict )
        self.assertTrue( len( data ) > 1 )

        expected = ['00232', '00983', '05204', '01100']
        result = self.gm.protein_maps['hsa:10']

        self.assertTrue( set(expected) == set(result) ) 


    def test_map_ecs( self ):

        data = self.gm._map_ecs()

        self.assertTrue( type( data ) is dict )
        self.assertTrue( len( data ) > 1 )

        expected = ['2.4.1.17', '1.2.1.31', '1.1.1.22', '1.2.1.47', '1.13.99.1', '1.2.1.8', '1.2.1.3', '3.1.1.17']
        result = self.gm.map_ecs['00053']

        self.assertTrue( set(expected) == set(result) )

    def test_organism_maps( self ):

        data = self.gm._organism_maps()

        self.assertTrue( type( data ) is list )
        self.assertTrue( len( data ) > 1 )

        self.assertEqual( len(self.gm._organism_maps()), 318 )


    def test_all_ec_numbers( self ):

        data = self.gm.all_ec_numbers

        self.assertTrue( type( data ) is list )
        self.assertTrue( len( data ) > 1 )

        self.assertEqual( len(data), 1153 )

    def test_all_incomplete_ec_numbers( self ):

        data = self.gm.all_incomplete_ec_numbers()

        incomplete = False

        for d in data:
            if re.findall('-', d ):
                incomplete = True
                break

        self.assertTrue( type( data ) is list )
        self.assertTrue( len( data ) > 1 )
        self.assertTrue( incomplete )


    def test_ec_by_protein_identification( self ):

        data = self.gm.ec_by_protein_identification( 'hsa:10327' )

        self.assertTrue( type( data ) is list )
        self.assertTrue( len( data ) > 0 )


if __name__ == "__main__":
    unittest.main()
