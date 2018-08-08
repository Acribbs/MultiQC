#!/usr/bin/env python

""" MultiQC module to parse the output of featurecounts for bespoke plotting"""

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
    """custom featurecounts output module, parses output logs. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name="Counting over genomic features", anchor="feature_counting",
        href='',
        info="This output plots the counts of reads to features, given your mapped bam file")

        # Init data
        self.feature_data = dict()
        
        # parse nm log
        for f in self.find_log_files('feature_counting/featuresmall'):
            self.feature_data[f['s_name'].replace(".feature_small","")] = self.parse_feature_logs(f['f'])

        if len(self.feature_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.feature_data)))

        # nm plot
        self.add_section(plot = self.feature_plot())

    def parse_feature_logs(self, f):
        data = {}

        n = 0

        for l in f.splitlines():
            n += 1
            s = l.split()
            if n <= 2:
                pass
            else:
                data[s[0]] = s[1]
        
        return data

    def feature_plot(self):
        """Make a barplot to show the proportion of found features """

        keys = list()

        # get the names of the features
        for s_name in self.feature_data:
            for features in self.feature_data[s_name]:
                keys.append(features)
        
        config = {
            'id': 'feature counts',
            'title': 'features'
            }
        
        return bargraph.plot(self.feature_data, keys, config)
