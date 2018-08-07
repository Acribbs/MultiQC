#!/usr/bin/env python

""" MultiQC module to parse the output of bam2stats-mapping quality of alignments"""

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ bam2stats module, parses output logs from mapq files. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name="bam2stats-mapping-quality", anchor="bam2stats",
        href='https://github.com/cgat-developers/cgat-apps/blob/master/CGAT/tools/bam2stats.py',
        info="bam2stats is a python script that will perform stats accross a bam file")

        # Init data
        self.bam2stats_data = dict()
        
        # parse nm log
        for f in self.find_log_files('bam2stats_mapq/mapq'):
            self.bam2stats_data[f['s_name'].replace(".readstats.mapq","")] = self.parse_bam2stats_logs(f['f'])

        # Filter to strip out samples
        self.bam2stats_data = self.ignore_samples(self.bam2stats_data)

        if len(self.bam2stats_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bam2stats_data)))

        # Write parsed report data to a file
        self.write_data_file(self.bam2stats_data, "multiqc_bam2stats_mapq")

        # Basic stats table
        self.bam2stats_general_stats_table()

        # nm plot
        self.add_section(plot = self.bam2stats_mapq_plot())



    def parse_bam2stats_logs(self, f):
        data = {}

        total = 0
        
        for l in f.splitlines():
            s = l.split()
            data[s[0]] = s[1]

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


    def bam2stats_mapq_plot(self):
        """Make the charts to show the alignments"""
        keys = list()
        pdata = OrderedDict()
        
        # get a list of the mapping quality names
        for s_name in self.bam2stats_data:
            for quality in self.bam2stats_data[s_name]:
                keys.append(quality)

        # now I have mapping quality names add the counts to it
        for s_name in self.bam2stats_data:
            pdata[s_name] = OrderedDict()
            for k in keys:
                if k == "mapq":
                    pass
                else:
                    pdata[s_name][k] = int(self.bam2stats_data[s_name][k])

        config = {
                'id': 'bam2stats-mapping-quality',
                'title': 'bam2stats mapping quality of alignments',
                'ylab': '# mapped alignments',
                'xlab': 'mapping quality score',
                'categories': True,
                'tt_label': '<strong>{point.category}:</strong> {point.y:.2f}'
            }
        
        return linegraph.plot(pdata, config)
