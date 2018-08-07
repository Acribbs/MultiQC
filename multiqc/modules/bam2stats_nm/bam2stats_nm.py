#!/usr/bin/env python

""" MultiQC module to parse the output of bam2stats-mismatches"""

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ bam2stats module, parses output logs. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name="bam2stats-mismatches", anchor="bam2stats",
        href='https://github.com/cgat-developers/cgat-apps/blob/master/CGAT/tools/bam2stats.py',
        info="bam2stats is a python script that will perform stats accross a bam file")

        # Init data
        self.bam2stats_data = dict()
        
        # parse nm log
        for f in self.find_log_files('bam2stats_nm/nm'):
            self.bam2stats_data[f['s_name'].replace(".readstats.nm","")] = self.parse_bam2stats_logs(f['f'])

        # Filter to strip out samples
        self.bam2stats_data = self.ignore_samples(self.bam2stats_data)

        if len(self.bam2stats_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bam2stats_data)))

        # Write parsed report data to a file
        self.write_data_file(self.bam2stats_data, "multiqc_bam2stats_nm")

        # Basic stats table
        self.bam2stats_general_stats_table()

        # nm plot
        self.add_section(plot = self.bam2stats_nm_plot())



    def parse_bam2stats_logs(self, f):
        data = {}

        total = 0
        
        for l in f.splitlines():
            s = l.split()
            if s[1] == "alignments":
                pass
            else:
                if int(s[0]) > 0:
                    total = int(total) + int(s[1])
                    data[s[0]] = s[1]
                else:
                    data[s[0]] = s[1]
        data['total_mismatches'] = total
        total = 0
        return data

    def bam2stats_general_stats_table(self):
        """
        This takes the parsed stats from kallisto report and then adds it to
        the basic stats table at the top of the report
        """
        headers = OrderedDict()

        headers['total_mismatches'] = {
            'title': 'mismatches',
            'description': 'The total number of mismatches',
            'scale': 'RdYlGn'
        }
        self.general_stats_addcols(self.bam2stats_data, headers)


    def bam2stats_nm_plot(self):
        """Make the charts to show the alignments"""
        keys = OrderedDict()
        keys['0'] = {'color': '#34495E', 'name': '0'}
        keys['1'] = {'color': '#FF5733', 'name': '1'}
        keys['2'] = {'color': '#FFC300', 'name': '2'}
        keys['3'] = {'color': '#DAF7A6', 'name': '3'}
        keys['4'] = {'color': '#900C3F', 'name': '4'}
        keys['5'] = {'color': '#581845', 'name': '5'}


        config = {
            'id': 'bam2stats missmatches',
            'title': 'mismatches'
            }
        
        return bargraph.plot(self.bam2stats_data, keys, config)
