

from Bio import SeqIO
import concurrent.futures
import time
import os
import shutil
import pickle
import gc
import random
from ete3 import Tree
import xml.etree.ElementTree as ET
import copy

# import networkx as nx
# import matplotlib.pyplot as plt
from . import _wrappers
from . import _utils_subhog
from . import _utils_frag_SO_detection
from ._hog_class import HOG

from ._wrappers import logger

low_so_detection = True # detection of proteins with low species overlap score in gene tree
fragment_detection = True  # this also need to be consistent in _hog_class.py
keep_subhog_each_pickle = False
inferhog_resume_subhog = True
inferhog_max_workers_num = 6 # how parallel to work in a species tree, (won't help near the root)
inferhog_min_hog_size_xml = 2     # by setting this as 1, pyham won't work on xml output.
orthoxml_v03 = True


def read_infer_xml_rhogs_batch(rhogid_batch_list, inferhog_concurrent_on, pickles_rhog_folder, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a list of rootHOGs
    """
    logger.debug("Inferring subHOGs for  "+str(len(rhogid_batch_list))+"rootHOGs started.")
    logger.debug("we are not reporting single tone hogs in the output xml. You may check this _config.inferhog_min_hog_size_xml.")
    hogs_rhog_xml_len_batch = []
    for rhogid in rhogid_batch_list:
        hogs_rhogs_xml_len = read_infer_xml_rhog(rhogid, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
        hogs_rhog_xml_len_batch.append(hogs_rhogs_xml_len)

    return hogs_rhog_xml_len_batch


def read_infer_xml_rhog(rhogid, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a  rootHOGs
    """

    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    if not os.path.exists(pickles_subhog_folder):
        os.makedirs(pickles_subhog_folder)

    # if (_config.gene_trees_write or _config.msa_write) and not os.path.exists("./genetrees"):
    #     os.makedirs("./genetrees")

    logger.debug("\n" + "==" * 10 + "\n Start working on root hog: " + rhogid + ". \n")
    rhog_i_prot_address = rhogs_fa_folder + "/HOG_" + rhogid+ ".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    logger.debug("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")
    # the file "species_tree_checked.nwk" is created by the check_input.py
    (species_tree) = _utils_subhog.read_species_tree(conf_infer_subhhogs.species_tree)

    (species_tree, species_names_rhog, prot_names_rhog) = _utils_subhog.prepare_species_tree(rhog_i, species_tree, rhogid)
    species_names_rhog = list(set(species_names_rhog))
    logger.debug("Number of unique species in rHOG " + rhogid + " is " + str(len(species_names_rhog)) + ".")

    if inferhog_concurrent_on:  # for big HOG we use parallelization at the level taxonomic level using concurrent
        hogs_a_rhog_num = infer_hogs_concurrent(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
    else:
        hogs_a_rhog_num = infer_hogs_for_rhog_levels_recursively(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
    # Output value hogs_a_rhog_num  is an integer= length. We save the output as pickle file at each taxonomic level.

    #####  Now read the final pickle file for this rootHOG
    root_node_name = species_tree.name
    pickle_subhog_file = pickles_subhog_folder + str(root_node_name) + ".pickle"
    with open(pickle_subhog_file, 'rb') as handle:
        hogs_a_rhog = pickle.load(handle)

    if not keep_subhog_each_pickle:
        shutil.rmtree(pickles_subhog_folder)

    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        if len(hog_i._members) >= inferhog_min_hog_size_xml:
            # could be improved   # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_a_rhog_xml_raw = hog_i.to_orthoxml()    # <generef  >      <paralg object >
            if orthoxml_v03 and 'paralogGroup' in str(hogs_a_rhog_xml_raw) :
                # in version v0.3 of orthoxml, there shouldn't be any paralogGroup at root level. Let's put them inside an orthogroup should be in
                hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(hog_i._hogid)})
                property_element = ET.SubElement(hog_elemnt, "property", attrib={"name": "TaxRange", "value": str(hog_i._tax_now)})
                hog_elemnt.append(hogs_a_rhog_xml_raw)
                hogs_a_rhog_xml = hog_elemnt
            else:
                hogs_a_rhog_xml = hogs_a_rhog_xml_raw
            hogs_rhogs_xml.append(hogs_a_rhog_xml)
        else:
            logger.debug("we are not reporting due to fastoma signleton hog |*|  " + str(list(hog_i._members)[0]))
            if len(hog_i._members)>1:
                logger.warning("issue 166312309 this is not a singleton "+str(hog_i._members))


    pickles_rhog_file = pickles_rhog_folder + '/file_' + rhogid + '.pickle'
    with open(pickles_rhog_file, 'wb') as handle:
        # dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger.debug("All subHOGs for the rootHOG as OrthoXML format is written in " + pickles_rhog_file)
    # to see orthoxml as string, you might need to do it for different idx
    # idx=0; from xml.dom import minidom; import xml.etree.ElementTree as ET; minidom.parseString(ET.tostring(hogs_rhogs_xml[idx])).toprettyxml(indent="   ")
    del hogs_a_rhog  # to be memory efficient
    gc.collect()
    hogs_rhogs_xml_len = len(hogs_rhogs_xml)
    return hogs_rhogs_xml_len


def infer_hogs_concurrent(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG using multi-threading (in parallel) on different taxanomic levels of species tree
    """

    pending_futures = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers= inferhog_max_workers_num) as executor:
        for node in species_tree.traverse(strategy="preorder"):
            node.dependencies_fulfilled = set()  # a set
            # node.infer_submitted = False
            if node.is_leaf():
                future_id = executor.submit(singletone_hog_, node, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
                # singletone_hog_(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
                # {<Future at 0x7f1b48d9afa0 state=finished raised TypeError>: 'KORCO_'}
                pending_futures[future_id] = node.name
                # node.infer_submitted = True
        # to keep, have a family
        while len(pending_futures) > 0:
            time.sleep(0.01)
            future_id_list = list(pending_futures.keys())
            for future_id in future_id_list:
                if future_id.done(): # done means finished, but may be unsecussful.

                    species_node_name = pending_futures[future_id]
                    logger.debug("checking for rootHOG id "+str(rhogid)+" future object is done for node " +str(species_node_name))
                    future_id.result()
                    #logger.debug("the result of future is  " + str(future_id.result()))

                    del pending_futures[future_id]
                    species_node = species_tree.search_nodes(name=species_node_name)[0]
                    parent_node = species_node.up
                    if not parent_node:  # we reach the root
                        # assert len(pending_futures) == 0, str(species_node_name)+" "+rhogid_
                        assert species_node.name == species_tree.name
                        break
                    parent_node.dependencies_fulfilled.add(species_node_name)  # a set

                    childrend_parent_nodes = set(node.name for node in parent_node.get_children())
                    if parent_node.dependencies_fulfilled == childrend_parent_nodes:
                        #  if not parent_node.infer_submitted:
                        future_id_parent = executor.submit(infer_hogs_this_level, parent_node, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs)
                        # parent_node.infer_submitted = True
                        # future_id_parent= parent_node.name+"aaa"
                        pending_futures[future_id_parent] = parent_node.name
                        # for future_id:  del pending_futures[future_id] i need another dictionary the other way arround to removes this futures

    return len(pending_futures) + 1


def infer_hogs_for_rhog_levels_recursively(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG using recursive function to traverse species tree (different taxanomic levels)
    """

    if sub_species_tree.is_leaf():
        singletone_hog_out = singletone_hog_(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
        #  out 1 =  succesful
        return singletone_hog_out
    children_nodes = sub_species_tree.children

    for node_species_tree_child in children_nodes:
        hogs_chrdn = infer_hogs_for_rhog_levels_recursively(node_species_tree_child, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
        # hogs_chrdn should be 1 hogs_chrdn_list.extend(hogs_chrdn)
    infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs)
    # hogs_this_level_list should be one
    return infer_hogs_this_level_out


def singletone_hog_(node_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    create subHOGs for leaves of tree (species level). Each protein is a subHOG.
    """

    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    this_level_node_name = node_species_name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    # logger.debug(" ** inferhog_resume_subhog is " + str(inferhog_resume_subhog))
    if inferhog_resume_subhog:

        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
        # open already calculated subhogs , but not completed till root in previous run
        if os.path.exists(pickle_subhog_file):
            if os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
                with open(pickle_subhog_file, 'rb') as handle:
                    # i don't even need to open this even
                    # is output of pickle.load(handle) is chlired or this level ?
                    # todo I think I don't need to read the pickle file, we could do it to check it has enough
                    hogs_this_level_list = pickle.load(handle) #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                    if hogs_this_level_list:
                        logger.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                        return len(hogs_this_level_list)
                    else:
                        logger.debug(" Issue  1238510: the pickle file for single tone is empty "+ str(hogs_this_level_list)+" " + rhogid)

    # logger.debug("reading protien / singletone HOG of  " + str(this_level_node_name))
    rhog_i_prot_address = rhogs_fa_folder +"/HOG_"+rhogid+".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    species_names_rhog_nonuniq = [seq.id.split("||")[1] for seq in rhog_i]
    prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog_nonuniq)) if
                                 species_names_rhog_nonuniq[idx] == node_species_name]
    rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]
    hogs_this_level_list = []
    for prot in rhog_part:
        hog_leaf = HOG(prot, node_species_name, rhogid)  # node_species_tree.name
        hogs_this_level_list.append(hog_leaf)
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)+".pickle"
    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs is written in pickle file.")

    return len(hogs_this_level_list)


def read_children_hogs(node_species_tree, rhogid, pickles_subhog_folder_all):
    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    # pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
    # if _config.inferhog_resume_subhog:  # open already calculated subhogs , but may not be completed in previous run
    #     if os.path.exists(pickle_subhog_file):
    #         if os.path.getsize(pickle_subhog_file) > 3:  # bytes
    #             with open(pickle_subhog_file, 'rb') as handle:
    #                 # i don't even need to open this even
    #                 hogs_children_level_list = pickle.load(handle)
    #                 if hogs_children_level_list:
    #                     return hogs_children_level_list

    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list = []
    for child_name in children_name:
        pickle_subhog_file = pickles_subhog_folder + str(child_name) + ".pickle"
        with open(pickle_subhog_file, 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))  # when there is an error somewhere else, probably with paralelization, for big rootHOG, the root reason of the error won't shown up. but you stope here at worst which is late.
            #  todo, check file exist, how to handle if not
    logger.debug("Finding hogs for rhogid: " + rhogid+ ", for taxonomic level:" + str(
        node_species_tree.name) + " for species sub-tree:\n  " + str(node_species_tree.write(format=1, format_root_node=True)) + "\n")
    return hogs_children_level_list




def infer_hogs_this_level(node_species_tree, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG at a taxanomic level
    """

    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    if node_species_tree.is_leaf():
        logger.warning("issue 1235,it seems that there is only on species in this tree, singleton hogs are treated elsewhere "+ rhogid)

    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"

    # TODO arrage resume with nextflow and also for when read single_tone pickles
    if inferhog_resume_subhog:
        if os.path.exists(pickle_subhog_file) and os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
            with open(pickle_subhog_file, 'rb') as handle:
                # todo : do I really need to read the pickle file
                hogs_this_level_list = pickle.load(handle)  #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                if hogs_this_level_list:
                    logger.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                    return len(hogs_this_level_list)
                else:
                    logger.debug(" Issue  1238510: the pickle file for single tone is empty " + str(hogs_this_level_list) + " " + rhogid)

    hogs_children_level_list = read_children_hogs(node_species_tree, rhogid, pickles_subhog_folder_all)

    if len(hogs_children_level_list) == 1:
        hogs_this_level_list = hogs_children_level_list
        with open(pickle_subhog_file, 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return len(hogs_children_level_list)


    genetree_msa_file_addr = "HOG_"+rhogid+"_"+str(node_species_tree.name) #+ ".nwk" # genetrees
    if len(genetree_msa_file_addr) > 245:
        # there is a limitation on length of file name. I want to  keep it consistent ,msa and gene tree names.
        rand_num = random.randint(1, 100000)
        genetree_msa_file_addr = genetree_msa_file_addr[:245] +"randomID_"+ str(rand_num) # +".nwk"
        logger.debug("genetree_msa_file_addr was long, now truncated " + genetree_msa_file_addr)
    genetree_msa_file_addr+="_itr0" # warning we expect this variable always ends with an integer (best 1 digit)
    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list if len(hog._msa) > 0]
    # sub_msa_list_lowerLevel_ready = [ii for ii in sub_msa_list_lowerLevel_ready_raw if len(ii) > 0]
    logger.debug("Merging "+str(len(sub_msa_list_lowerLevel_ready))+" MSAs for rhog:"+rhogid+", level:"+str(node_species_tree.name))
    msa_filt_row_col = []
    prot_dubious_msa_list = []
    seq_dubious_msa_list = []
    merged_msa = ""
    if sub_msa_list_lowerLevel_ready:
        if len(sub_msa_list_lowerLevel_ready) > 1:
            merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready, genetree_msa_file_addr, conf_infer_subhhogs)
            if fragment_detection:
                prot_dubious_msa_list, seq_dubious_msa_list = _utils_frag_SO_detection.find_prot_dubious_msa(merged_msa, conf_infer_subhhogs)
        else:
            merged_msa = sub_msa_list_lowerLevel_ready   #  when only on  child, the rest msa is empty.
        logger.debug("All sub-hogs are merged, merged_msa "+str(len(merged_msa))+" "+str(len(merged_msa[0]))+" for rhog: "+rhogid+", taxonomic level:"+str(node_species_tree.name))
        (msa_filt_row_col, msa_filt_col, hogs_children_level_list, genetree_msa_file_addr) = _utils_subhog.filter_msa(merged_msa, genetree_msa_file_addr, hogs_children_level_list, conf_infer_subhhogs)
        # msa_filt_col is used for parent level of HOG. msa_filt_row_col is used for gene tree inference.
    else:
        logger.info("Issue 1455, merged_msa is empty " + rhogid + ", for taxonomic level:" + str(node_species_tree.name))

    if len(msa_filt_row_col) > 1 and len(msa_filt_row_col[0]) > 1:

        gene_tree_raw = _wrappers.infer_gene_tree(msa_filt_row_col, genetree_msa_file_addr, conf_infer_subhhogs)
        gene_tree=""
        try:
            gene_tree = Tree(gene_tree_raw + ";", format=0)   #
        except:
            try:
                gene_tree = Tree(gene_tree_raw + ";", format=0, quoted_node_names=True)  #
            except:
                logger.debug(" issue 22233198233 error"+str(gene_tree))
        logger.debug("Gene tree is inferred len "+str(len(gene_tree))+" rhog:"+rhogid+", level: "+str(node_species_tree.name))

        if fragment_detection and len(gene_tree) > 2 and prot_dubious_msa_list:
            (gene_tree, hogs_children_level_list, merged_msa_new) = _utils_frag_SO_detection.handle_fragment_msa(prot_dubious_msa_list, seq_dubious_msa_list, gene_tree, node_species_tree, genetree_msa_file_addr, hogs_children_level_list, merged_msa, conf_infer_subhhogs)
        else:
            merged_msa_new = merged_msa

        # when the prot dubious is removed during trimming
        if len(gene_tree) > 1: # e.g. "('sp|O67547|SUCD_AQUAE||AQUAE||1002000005|_|sub10001':0.329917,'tr|O84829|O84829_CHLTR||CHLTR||1001000005|_|sub10002':0.329917);"
            (gene_tree, all_species_dubious_sd_dic, genetree_msa_file_addr) = _utils_subhog.genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr,conf_infer_subhhogs, hogs_children_level_list)

            if low_so_detection and all_species_dubious_sd_dic:
                (gene_tree, hogs_children_level_list, genetree_msa_file_addr) = _utils_frag_SO_detection.handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic, hogs_children_level_list, conf_infer_subhhogs)

            logger.debug("Merging sub-hogs for rhogid:"+rhogid+", level:"+str(node_species_tree.name))
            # the last element should be merged_msa not the trimmed msa, as we create new hog based on this msa
            hogs_this_level_list = merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid, merged_msa_new, conf_infer_subhhogs)
            # for i in hogs_this_level_list: print(i.get_members())
            logger.debug("After merging subhogs of childrens, "+str(len(hogs_this_level_list))+" subhogs are found for rhogid: "+rhogid+", for taxonomic level:"+str(this_level_node_name))

        else:
            hogs_this_level_list = hogs_children_level_list
    else:
        if msa_filt_row_col:
            logger.debug("warning id 13805: hogs_this_level_list is empty. msa_filt_row_col:"+str(len(msa_filt_row_col))+"*"+str(len(msa_filt_row_col[0]))+" !!")
        else:
            logger.debug("warning id 13806: msa_filt_row_col is empty." + str(len(msa_filt_row_col)) +"! ")

        hogs_this_level_list = hogs_children_level_list

    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # for hog in hogs_this_level_list:
    #     msa_rec_ids = [i.id for i in hog._msa]
    #     if set(hog._members) != set(msa_rec_ids):
    #         logger_hog.warning("issue 123601, the members are not matching with msa"+str(set(hog._members))+"  "+str(msa_rec_ids))
    # this could be becuase of sub-sampling


    return len(hogs_this_level_list)


def remove_protlist_hog_hierarchy_toleaves(hog_ii, protlist_to_remove):
    # import copy; a= copy.deepcopy(hog_ii)
    for subhog in hog_ii._subhogs:
        intersection_prots = list(subhog._members.intersection(set(protlist_to_remove)))
        if  intersection_prots:
            result_removing0 = remove_protlist_hog_hierarchy_toleaves(subhog, intersection_prots)

    result_removing = hog_ii.remove_protlist_from_hog(intersection_prots)
    subhogs_old = hog_ii._subhogs
    hog_ii._subhogs = [i for i in subhogs_old if len(i._members)]  # for getting rid of  empty subhogs #  hogID=HOG_0026884_sub12868,len=0, tax_least=Ardeidae

    return result_removing


def keep_prots_hog_hierarchy_toleaves(hog_ii, prots_tokeep):

    for subhog in hog_ii._subhogs:
        prots_intser= subhog._members.intersection((set(prots_tokeep)))
        if prots_intser:
            result_removing0 = keep_prots_hog_hierarchy_toleaves(subhog, prots_intser)
        else:
            subhog._members = set([])

    result_pruning = hog_ii.prune(prots_tokeep)
    subhogs_old = hog_ii._subhogs
    hog_ii._subhogs = [i for i in subhogs_old if len(i._members)]  # for getting rid of  empty subhogs #  hogID=HOG_0026884_sub12868,len=0, tax_least=Ardeidae

    return result_pruning


def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid, merged_msa, conf_infer_subhhogs):
    """
    merge subhogs based on the gene tree speciation node of gene tree by creating inter-HOG graph (implicitely )
    """
    hogs_children_level_updated = copy.deepcopy(hogs_children_level_list)
    prots_genetree_raw0 = [i.name for i in gene_tree.get_leaves()]
    prots_genetree_raw1 = [i.split("|_|")[0] for i in prots_genetree_raw0]
    if prots_genetree_raw1[0].startswith("'"):
        prots_genetree = [i[1:] for i in prots_genetree_raw1]  # since is split , there is no ' at the end. so -1] is not needed. # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005", gene tree is quoted, there are both ' and " !
    else:
        prots_genetree = prots_genetree_raw1
        # todo how about merged split genes, they have _|_ inside   (also for below leaves_name_raw0)
        # leaves_name =[] for name_i in leaves_name): leaves_name += name_i.split("_|_")

    hogs_this_level_list = []
    subHOG_ids_merged_done = []
    for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "processed") and n.processed==True):

        if not node.is_leaf() and node.name[0] == "S":
            leaves_name_raw0 = [i.name for i in node.get_leaves()]  # leaves names  with subhog id  'HALSEN_R15425||HALSEN|_|1352015793_sub10149',
            leaves_name_raw1 = [i.split("|_|")[0] for i in leaves_name_raw0]
            if leaves_name_raw1[0].startswith("'"):
                S_leaves_name = [i[1:] for i in leaves_name_raw1] # since is split, there is no ' at the end. so -1] is not needed.  # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005", gene tree is quoted, there are both ' and " !
            else:
                S_leaves_name = leaves_name_raw1

            subHOGs_to_merge = []
            for node_leave_name in S_leaves_name:
                for subHOG in hogs_children_level_updated:
                    subHOG_members = subHOG._members
                    if node_leave_name in subHOG_members:  #we merge a hog only once, unless it is split (the rest of hog will be merged in another speciation event
                        if subHOG not in subHOGs_to_merge and subHOG._hogid not in subHOG_ids_merged_done:
                            subHOGs_to_merge.append(subHOG)

            if len(subHOGs_to_merge)<2:
                continue

            for subhog in subHOGs_to_merge:
                subhog_members_intree = [i for i in subhog._members if i in prots_genetree]
                subhog_members_intree_notInSpeciaion = [i for i in subhog_members_intree if i not in S_leaves_name]
                subhog_members_intree_InSpeciaion = [i for i in subhog_members_intree if i in S_leaves_name]

                if subhog_members_intree_notInSpeciaion:  # a few of the proteins in the subhog (exist in gene tree) are in other part of the tree, due to dup probably
                    print("Inconsistency compared previous level, we'll split the hog"+str(subhog))  # otherwise all prots in this hog were in the same child subhog  in the previous level ->  subhog_rest =  subhog and subhog=empty

                    subhog_members_intree_notInSpeciaion_extended= subhog_members_intree_notInSpeciaion
                    # try to include unsampled prots (not in the tree), challange those that are removed due to low species overlap might come back
                    # to do so we use the subhog structure of previous level
                    # the followings are the children of children (children of subhogs)
                    subhog_childids_intree_notInSpeciaion = []
                    for prot in subhog_members_intree_notInSpeciaion:
                        for child_subhog in subhog._subhogs:
                            if prot in child_subhog._members and child_subhog not in subhog_childids_intree_notInSpeciaion:
                                subhog_childids_intree_notInSpeciaion.append(child_subhog)

                    prots_childids = [i for subhog in subhog_childids_intree_notInSpeciaion for i in subhog._members]
                    # such child subhog might have some proteins of the speciation node and some outside of specietion node (contradactory info, we do our best based on the latest information we have)
                    subhog_members_intree_notInSpeciaion_extended = [i for i in prots_childids if i not in subhog_members_intree_InSpeciaion]
                    subhog_members_rest = [i for i in subhog._members if i in subhog_members_intree_notInSpeciaion_extended]

                    hogs_children_level_updated.remove(subhog)
                    subhog_rest = copy.deepcopy(subhog)
                    subhog_hogid = subhog._hogid
                    subhog._hogid =    subhog_hogid+"|0"
                    subhog_rest._hogid=subhog_hogid+"|1"

                    for prot_ii in subhog_members_intree_notInSpeciaion_extended:
                        result_removing = _utils_frag_SO_detection.remove_prot_hog_hierarchy_toleaves(subhog, prot_ii)
                        if result_removing == 0:   # the hog is empty
                            print("warning result_removing is empty!"+str(subhog))
                    print("The subhog is shrunk to" + str(subhog))
                    results_keep = keep_prots_hog_hierarchy_toleaves(subhog_rest, subhog_members_rest)
                    print("The subhog_rest is here " + str(subhog_rest))
                    hogs_children_level_updated.append(subhog_rest)
                    hogs_children_level_updated.append(subhog)

                    # if set(prot_list_notintheSpeciaionNode_extended) == set(subhog._members):
                    #    logger.warning("issue 123550971 we are merging two subhogs that there was a dup, we couldn't split because proteins in the species node and outside of it were in the same child subhog (unsampled)   ")

            # add a check when subhog_rest is the last item in the gene tree, to include non-sampled genes
            taxnomic_range = node_species_tree.name
            num_species_tax_speciestree = len(node_species_tree.get_leaves()) # num_species_tax   is the number of species exist in the species tree at this clade
            subHOG_ids_merged_done += [hg._hogid for hg in subHOGs_to_merge]
            HOG_this_node = HOG(set(subHOGs_to_merge), taxnomic_range, rhogid, merged_msa, num_species_tax_speciestree, conf_infer_subhhogs)
            if len(HOG_this_node._msa) == 1:
                logger.warning("issue 1258313"+str(HOG_this_node)+str(HOG_this_node._msa)+" "+node.name  )
            hogs_this_level_list.append(HOG_this_node)


                #  I don't need to traverse deeper in this clade
            node.processed = True  # print("?*?*  ", node.name)

    subhogs_id_children_assigned = [child_subhog._hogid for hog in hogs_this_level_list for child_subhog in hog._subhogs]
    for subHOG in hogs_children_level_updated:  # for the single branch  ( D include subhog and S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    for hog_j in hogs_this_level_list:
        hog_j._tax_now = node_species_tree.name

    return hogs_this_level_list
