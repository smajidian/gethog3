#!/usr/bin/env python3
# from pyoma project

from __future__ import unicode_literals, division
from builtins import str
import lxml.etree as etree
import os
import errno
import logging
import re
import argparse
import sys
logger = logging.getLogger(__name__)


## the only change compared to pyoma version is to remove int constratin on HOG id .


class OrthoXMLSplitter(object):
    """Convert orthoxml files with several families.

    This class provides the means to extract a subset of root HOGs (i.e.
    families) into a new output orthoxml file, or to split it and create
    for each family an individual file.

    The object should be instantiated with the input orthoxml file and
    optionally a cache_dir argument where the output orthoxml files will
    be stored. This later parameter can be overwritten in the __call__
    method call that does the work.

    .. note::

       Calls to the splitter will remove the created families from the
       loaded input file, so subsequent calls that contain a family in
       common will miss them from the second call onwards.


    :Example:

      splitter = OrthoXMLSplitter("data.orthoxml", cache_dir="./splits")
      splitter()

    will create files HOGxxxxxx.orthoxml in the ./splits directory."""

    def __init__(self, xml_file, cache_dir=None, release_char=None):
        self.xml_file = xml_file
        if cache_dir is not None:
            self._assert_cache_dir(cache_dir)
        if release_char is None:
            self.release_char = ""
        elif re.match(r"^[A-Z]?$", release_char):
            self.release_char = release_char
        else:
            raise ValueError(
                "unexpected value for release_char: '{}'. Needs to be a single capital ascii letter".format(
                    release_char
                )
            )
        logger.info("loading xml file {}...".format(xml_file))
        parser = etree.XMLParser(remove_blank_text=True)
        self.Etree_XML = etree.parse(self.xml_file, parser=parser)
        self.Etree_root = self.Etree_XML.getroot()
        logger.info("building lookup table for genes")
        self.gene_lookup = {gene.get("id"): gene for gene in self._iter_gene_elements()}
        logger.info("init of OrthoXMLSplitter finished")

    def _assert_cache_dir(self, cache_dir):
        # Ensure existance of cache directory (py2 compat)
        try:
            os.makedirs(cache_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(cache_dir):
                pass
            else:
                raise
        self.cache_dir = cache_dir

    def _iter_gene_elements(self):
        """This method is a faster version of xpath '//ns:gene'.

        It iterates the element in sequential order"""
        for node in self.Etree_root:
            if node.tag == "{http://orthoXML.org/2011/}species":
                for gene in node.iter("{http://orthoXML.org/2011/}gene"):
                    yield gene

    def _iter_toplevel_groups(self):
        """This method yields all the root hogs sequentially."""
        for node in self.Etree_root:
            if node.tag == "{http://orthoXML.org/2011/}groups":
                for root_hog in node:
                    yield root_hog

    def __call__(
        self,
        hogs_to_extract=None,
        single_hog_files=False,
        basename=None,
        cache_dir=None,
    ):
        """Split/extract hogs from orthoxml file based on root hogs ids.

        Split the input orthoxml or extract a subset of root hogs. If no
        argument is passed, one orthoxml file per root hog is created,
        named as 'HOGxxxxxx.orthoxml', where xxxxxx is the numeric id of
        each hog.

        The set of root hogs to be extracted can be limited by specifying
        a subset of hog ids in the hogs_to_extract parameter. If
        single_hog_files is set to true, each of these hogs will be converted
        into a single orthoxml file named as explained above. If single_hog_files
        is set to false, the whole subset of hogs will be stored in one
        orthoxml file named as specified in `basename`.

        The file(s) will be stored in the cache_dir folder which can be
        specified in the constructor or overwritten as an argument in
        this method.

        :param hogs_to_extract: list or set that contains the set of root
            hogs to be extracted. If set to None, all hogs are extracted.
        :param bool single_hog_files: whether or not to build one orthoxml
            file for all the selected hogs or individual ones.
        :param str basename: name of the output file if a subset of hogs
            is extracted into a single file.
        :param str cache_dir: folder where to store the output files.
        """
        if cache_dir is not None:
            self._assert_cache_dir(cache_dir)
        elif self.cache_dir is None:
            raise RuntimeError("cache dir to output files to is not set")

        if single_hog_files:
            if hogs_to_extract is None:
                raise RuntimeError(
                    "useless to extract all hogs into single output file"
                )
            if basename is None or not isinstance(basename, (str, bytes)):
                raise ValueError("basename needs to be specified: {}".format(basename))
            ogs = [
                og
                for og in self._iter_toplevel_groups()
                if og.get("id") in hogs_to_extract
            ]
            fn = os.path.join(self.cache_dir, basename)
            logger.info("extracting {:d} hogs into {:s}".format(len(ogs), fn))
            self.create_new_orthoxml(fn, ogs)
        else:
            for og in self._iter_toplevel_groups():
                if hogs_to_extract is None or og.get("id") in hogs_to_extract:
                    hog_nr = og.get("id")
                    hog_id = hog_nr+".orthoxml" #"HOG{:07d}.orthoxml".format(hog_nr)
                    fname = os.path.join(self.cache_dir, hog_id)
                    logger.info("extracting {} into {}".format(hog_id, fname))
                    self.create_new_orthoxml(fname, [og])

    def iter_generefs_in_og(self, og_node):
        for node in og_node.iterdescendants("{http://orthoXML.org/2011/}geneRef"):
            yield node

    def get_gene_via_generef(self, genesref_ids):
        genesref_ids = set(genesref_ids)
        return [self.gene_lookup[gene_id] for gene_id in genesref_ids]

    def create_new_orthoxml(self, fn, OGs):
        """create a new orthoxml file for the passed orthologGroup elements.

        :param fn: the filename of the output file. The path needs to exists
            prior to calling this method.
        :param OGs: the orthologGroup elements that should be included in the
            new output file."""
        # Get element to store
        for og_node in OGs:
            gene_ids = [
                gene_ref_elem.get("id")
                for gene_ref_elem in self.iter_generefs_in_og(og_node)
            ]
        gene_els = self.get_gene_via_generef(gene_ids)

        # Get all information to store
        zoo = {}  # <- {key:sp_etree || value: {key:db_el || values:[list_genes]}}
        for gene_el in gene_els:  # <- for all gene el
            db_el = gene_el.getparent().getparent()
            sp_el = db_el.getparent()
            if sp_el in zoo.keys():  # <- if species already visited
                if db_el in zoo[sp_el].keys():  # <- if db already visited so add gene
                    zoo[sp_el][db_el].append(gene_el)
                else:  # <- if db not visited so add db,genes
                    zoo[sp_el][db_el] = []
                    zoo[sp_el][db_el].append(gene_el)
            else:  # <- if species not visited so add sp,db,gene
                zoo[sp_el] = {}
                zoo[sp_el][db_el] = []
                zoo[sp_el][db_el].append(gene_el)

        etree_2_dump = etree.Element("orthoXML", nsmap=self.Etree_root.nsmap)
        for attr, value in self.Etree_root.items():
            etree_2_dump.set(attr, value)

        for species_el in zoo.keys():
            species_xml = etree.Element("species")
            for attr, value in species_el.items():
                species_xml.set(attr, value)
            etree_2_dump.insert(0, species_xml)

            for db_el in zoo[species_el].keys():
                # Add <database> into <species>
                database_xml = etree.SubElement(species_xml, "database")
                for attr, value in db_el.items():
                    database_xml.set(attr, value)

                # Add <genes> TAG into <database>
                genes_xml = etree.SubElement(database_xml, "genes")

                # Fill <genes> with <gene>
                for gene_el in zoo[species_el][db_el]:
                    gene_xml = etree.SubElement(genes_xml, "gene")
                    for attr, value in gene_el.attrib.items():
                        gene_xml.set(attr, value)

        groupsxml = etree.SubElement(etree_2_dump, "groups")
        for og_et in OGs:
            if not og_et.get("id").startswith("HOG:{:s}".format(self.release_char)):
                og_et.set(
                    "id",
                    og_et.get("id"),
                    #"HOG:{:s}{:07d}".format(self.release_char, int(og_et.get("id"))),
                )
            groupsxml.append(og_et)

        tree = etree.ElementTree(etree_2_dump)
        tree.write(
            fn, xml_declaration=True, encoding="utf-8", method="xml", pretty_print=True
        )
        
### main
def main(argv):
    parser = argparse.ArgumentParser(description="Split OrthoXML file into HOGs")
    parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="OrthoXML file")
    parser.add_argument("-o", "--outFile", dest="outFile", required=False, default='./splits/', help='split OrthoXML destination path')
    args = parser.parse_args()

    splitter = OrthoXMLSplitter(args.inFile, cache_dir=args.outFile)
    splitter()

if __name__ == '__main__':
    main(sys.argv)

