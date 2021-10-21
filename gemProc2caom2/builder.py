# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  : 4 $
#
# ***********************************************************************
#

import logging
from os.path import basename
from urllib.parse import urlparse

from cadctap import CadcTapClient
from caom2utils import fits2caom2
from caom2pipe import client_composable as clc
from caom2pipe import manage_composable as mc
from caom2pipe import name_builder_composable as nbc
from gem2caom2 import gem_name
from gem2caom2 import external_metadata as em
from gem2caom2.obs_file_relationship import repair_data_label


__all__ = ['CADC_SCHEME', 'COLLECTION', 'GemProcBuilder', 'GemProcName']


COLLECTION = 'GEMINICADC'
CADC_SCHEME = 'cadc'


class GemProcName(mc.StorageName):
    def __init__(self, entry):
        if entry.startswith('vos'):
            self._vos_uri = entry
            self._file_name = basename(urlparse(self._vos_uri).path)
            self._file_id = gem_name.GemName.remove_extensions(
                self._file_name
            )
        else:
            self._file_name = basename(entry)
            self._file_id = gem_name.GemName.remove_extensions(
                self._file_name
            )
        super(GemProcName, self).__init__(
            fname_on_disk=self._file_name,
            collection=COLLECTION,
            compression='',
            scheme=CADC_SCHEME,
            entry=entry,
        )
        self._source_names = [entry]
        self._destination_uris = [self.file_uri]
        self.fname_on_disk = self._file_name
        self._product_id = self._file_id
        # must use the builder
        self._obs_id = None
        self._logger = logging.getLogger(__name__)

    @property
    def file_id(self):
        return self._file_id

    @property
    def file_name(self):
        return self._file_name

    @property
    def prev(self):
        return '{}.jpg'.format(self._file_id)

    @property
    def product_id(self):
        return self._product_id

    @property
    def thumb(self):
        return '{}_th.jpg'.format(self._file_id)

    def is_valid(self):
        # over-ride self._obs_id dependency
        return True


class GemProcBuilder(nbc.StorageNameBuilder):
    def __init__(self, config):
        super().__init__()
        self._config = config
        self._subject = clc.define_subject(self._config)
        self._sc2_client = CadcTapClient(
            subject=self._subject, resource_id='ivo://cadc.nrc.ca/sc2tap'
        )
        self._prod_client = CadcTapClient(
            subject=self._subject, resource_id='ivo://cadc.nrc.ca/ams/gemini'
        )
        self._connected = mc.TaskType.SCRAPE not in config.task_types
        self._logger = logging.getLogger(self.__class__.__name__)

    def build(self, entry):
        self._logger.debug(f'Begin build for {entry}')
        temp = urlparse(entry)
        file_name = basename(temp.path)
        gem_proc_name = GemProcName(entry=entry)
        gem_proc_name.obs_id = self._get_obs_id(temp, file_name, entry)
        self._logger.debug(f'End build with {gem_proc_name}')
        return gem_proc_name

    def _get_obs_id(self, temp, file_name, entry):
        """
        These files are not available from archive.gemini.edu, so
        only ask for their metadata from CADC.
        """
        self._logger.debug(f'Begin _get_obs_id for file_name {file_name}')
        metadata = None
        if self._connected:
            if self._config.use_local_files:
                self._logger.debug(f'Check local {file_name}')
                metadata = em.defining_metadata_finder._check_local(file_name)
            if metadata is None and temp is not None and temp.scheme == 'vos':
                self._logger.debug('Check vos')
                metadata = self._get_obs_id_from_vos(entry)
            if metadata is None:
                # why the old collection name? Because it's better to
                # retrieve the metadata from the old sc2 collection than
                # by retrieving a header, and beat up CADC instead of
                # archive.gemini.edu
                original_client = em.defining_metadata_finder._tap_client
                try:
                    self._logger.debug(f'Check caom2 collection {COLLECTION}')
                    # uri = mc.build_uri(COLLECTION, file_name, CADC_SCHEME)
                    em.defining_metadata_finder._tap_client = (
                        self._prod_client
                    )
                    for uri in [
                        f'gemini:GEM/{file_name}',
                        f'gemini:GEMINI/{file_name}',
                    ]:
                        metadata = em.defining_metadata_finder._check_caom2(
                            uri, COLLECTION
                        )
                        if metadata is not None:
                            break
                    if metadata is None:
                        self._logger.debug(
                            f'Check caom2 collection GEMINIPROC'
                        )
                        # uri = mc.build_uri('GEMINI', file_name)
                        uri = f'ad:GEMINI/{file_name}'
                        em.defining_metadata_finder._tap_client = (
                            self._sc2_client
                        )
                        metadata = em.defining_metadata_finder._check_caom2(
                            uri, 'GEMINIPROC'
                        )
                finally:
                    em.defining_metadata_finder._tap_client = original_client
        else:
            self._logger.debug('Check unconnected local')
            metadata = em.defining_metadata_finder._check_local(file_name)
        if metadata is None:
            raise mc.CadcException(f'No metadata for {file_name}')
        if metadata.data_label is not None:
            metadata.data_label = repair_data_label(
                file_name, metadata.data_label
            )
        self._logger.debug(f'End _get_obs_id')
        return metadata.data_label

    def _get_obs_id_from_vos(self, entry):
        self._logger.debug(f'Begin get_obs_id_from_vos for {entry}.')
        headers = fits2caom2.get_vos_headers(entry, self._subject)
        if headers is None or len(headers) == 0:
            raise mc.CadcException(
                f'Could not get metadata from {self._vos_uri}'
            )
        else:
            data_label = mc.get_keyword(headers, 'DATALAB')
            instrument = mc.get_keyword(headers, 'INSTRUME')
            if data_label is None or instrument is None:
                raise mc.CadcException(
                    f'Could not find DATALAB {data_label} or INSTRUME '
                    f'{instrument} in {entry}.'
                )
            metadata = em.DefiningMetadata(instrument, data_label)
        self._logger.debug('End get_obs_id_from_vos.')
        return metadata
