# -*- coding: utf-8 -*-
"""
Created on Wed May 24 22:41:44 2023

@author: dkane
"""

import os
import pandas as pd
import struct

from tqdm import tqdm
from tkinter import filedialog
from Semi_ATE.STDF import utils
from Semi_ATE.STDF import (ATR, BPS, CDR, CNR, DTR, EPS, FAR, FTR, GDR, HBR, MIR, MPR, MRR, 
                NMR, PCR, PGR, PIR, PLR, PMR, PRR, PSR, PTR, RDR, SBR, SDR, SSR, 
                STR, TSR, VUR, WIR, WCR, WRR)


'''
Write STDF flow:
    1) read records from file into index dict
    2) make changes to index dict
    3) write records to stdf file
'''
class STDFFile:
    '''
        initialize 
    '''
    def __init__(self, fp, progress = False):
        self.fp = fp
        self.progress = progress
        self.index = self._index_stdf(self.fp)
        self.test_num_name = self._get_test_num_nam_dict()
        
        
    def _index_stdf(self, fp):
        index = {}
        offset = 0
        
        if utils.is_STDF(fp):
            endian, version = utils.endian_and_version_from_file(fp)
            index['version'] = version
            index['endian'] = endian
            index['records'] = {}
            index['indexes'] = {}
            index['parts'] = {}
            PIP = {} # parts in process
            PN = 1
            
            TS2ID = utils.ts_to_id(version)
            
            if self.progress:
                description = "Indexing STDF file '%s'" % os.path.split(fp)[1]
                index_progress = tqdm(total=os.path.getsize(fp), ascii=True, disable=not self.progress, desc=description, leave=False, unit='b')

            for _, REC_TYP, REC_SUB, REC in utils.check_records_from_file(fp):
                REC_ID = TS2ID[(REC_TYP, REC_SUB)]
                REC_LEN = len(REC)
                if REC_ID not in index['records']: index['records'][REC_ID] = []
                index['indexes'][offset] = REC
                index['records'][REC_ID].append(offset)
                if REC_ID in ['PIR', 'PRR', 'PTR', 'FTR', 'MPR']:
                    if REC_ID == 'PIR':
                        pir_HEAD_NUM, pir_SITE_NUM = self._get_head_and_site_num("PIR", REC)
                        if (pir_HEAD_NUM, pir_SITE_NUM) in PIP:
                            raise Exception("One should not be able to reach this point !")
                        PIP[(pir_HEAD_NUM, pir_SITE_NUM)] = PN
                        index['parts'][PN]=[]
                        index['parts'][PN].append(offset)
                        PN+=1
                    elif REC_ID == 'PRR':
                        prr_HEAD_NUM, prr_SITE_NUM = self._get_head_and_site_num("PRR", REC)
                        if (prr_HEAD_NUM, prr_SITE_NUM) not in PIP:
                            raise Exception("One should not be able to reach this point!")
                        pn = PIP[(prr_HEAD_NUM, prr_SITE_NUM)]
                        index['parts'][pn].append(offset)
                        del PIP[(prr_HEAD_NUM, prr_SITE_NUM)]
                    elif REC_ID == 'PTR':
                        ptr_HEAD_NUM, ptr_SITE_NUM = self._get_head_and_site_num("PTR", REC)
                        if (ptr_HEAD_NUM, ptr_SITE_NUM) not in PIP:
                            raise Exception("One should not be able to reach this point!")
                        pn = PIP[(ptr_HEAD_NUM, ptr_SITE_NUM)]
                        index['parts'][pn].append(offset)
                    elif REC_ID == 'FTR':
                        ftr_HEAD_NUM, ftr_SITE_NUM = self._get_head_and_site_num("FTR", REC)
                        if (ftr_HEAD_NUM, ftr_SITE_NUM) not in PIP:
                            raise Exception("One should not be able to reach this point!")
                        pn = PIP[(ftr_HEAD_NUM, ftr_SITE_NUM)]
                        index['parts'][pn].append(offset)
                    elif REC_ID == 'MPR':
                        mpr_HEAD_NUM, mpr_SITE_NUM = self._get_head_and_site_num("MPR", REC)
                        if (mpr_HEAD_NUM, mpr_SITE_NUM) not in PIP:
                            raise Exception("One should not be able to reach this point!")
                        pn = PIP[(mpr_HEAD_NUM, mpr_SITE_NUM)]
                        index['parts'][pn].append(offset)
                    else:
                        raise Exception("One should not be able to reach this point! (%s)" % REC_ID)
                if self.progress: index_progress.update(REC_LEN)
                offset += REC_LEN
        return index
                
    def _get_test_num_nam_dict(self):
        TEST_NUM_NAM = {}

        for tsr_offset in self.index['records']['TSR']:
            tsr = TSR(self.index['version'], self.index['endian'], self.index['indexes'][tsr_offset])
            TEST_NUM = tsr.get_value('TEST_NUM')
            TEST_NAM = tsr.get_value('TEST_NAM')
            TEST_TYP = tsr.get_value('TEST_TYP').upper()
            if TEST_NUM not in TEST_NUM_NAM:
                TEST_NUM_NAM[TEST_NUM] = []
            if (TEST_NAM, TEST_TYP) not in TEST_NUM_NAM[TEST_NUM]:
                TEST_NUM_NAM[TEST_NUM].append((TEST_NAM, TEST_TYP))
        return TEST_NUM_NAM
    
    # fast method to return head num and site num from raw_bytes
    def _get_head_and_site_num(self, rec_id, rec):
        if rec_id in ["PRR", "PIR"]:
            head_num = struct.unpack('b', rec[4:5])[0]
            site_num = struct.unpack('b', rec[5:6])[0]
        elif rec_id in ["FTR", "PTR", "MPR"]:
            head_num = struct.unpack('b', rec[8:9])[0]
            site_num = struct.unpack('b', rec[9:10])[0]
        else:
            raise Exception(f"method can't handle rec_id: {rec_id}")
        return head_num, site_num
    
    def write_stdf(self, fp, overwrite=True):
        if not overwrite:
            assert not os.path.isfile(fp), f"File already exists: {fp}"
        with open(fp, "wb") as stdf_file:
            for rec in self.index['indexes'].values():
                stdf_file.write(rec)
        
        
if __name__ == "__main__":
    # fp = r'C:/Users/dkane/OneDrive - Presto Engineering/Documents/python_scripts/semi ate stdf processing/stdf file/test files/123456_25_9_25_9__20230531_152904.stdf'
    # fp = r'C:/Users/dkane/OneDrive - Presto Engineering/Documents/python_scripts/semi ate stdf processing/stdf file/test files/5AIX5202-P102.std'
    fp = r'C:/Users/dkane/OneDrive - Presto Engineering/Documents/python_scripts/semi ate stdf processing/stdf file/test files/G4_MZMD15163_N19347.1_01_20230525.stdf'
    # fp = r'C:/Users/dkane/OneDrive - Presto Engineering/Documents/python_scripts/semi ate stdf processing/stdf file/test files/G4_TIA25133_N19347.1_01_20230530.stdf'
    
    # get indexed stdf file
    stdf = STDFFile(fp, progress = True)
    
    # build new stdf file name
    fp_no_ext = os.path.splitext(fp)[0]
    new_fp = fp_no_ext + "_EDITED.stdf"
    
    # write new stdf file
    stdf.write_stdf(new_fp, overwrite = True)
    
        
