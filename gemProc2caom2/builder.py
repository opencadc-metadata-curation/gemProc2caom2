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

from logging import getLogger
from os.path import basename
from urllib.parse import urlparse

from caom2pipe.manage_composable import (
    CadcException,
    get_keyword,
    StorageName,
)
from caom2pipe import name_builder_composable as nbc
from gem2caom2.obs_file_relationship import repair_data_label


__all__ = ['CADC_SCHEME', 'COLLECTION', 'GemProcBuilder', 'GemProcName']


COLLECTION = 'GEMINICADC'
CADC_SCHEME = 'cadc'


class GemProcName(StorageName):
    def __init__(self, entry):
        if entry.startswith('vos'):
            self._vos_uri = entry
            file_name = basename(urlparse(self._vos_uri).path)
        else:
            file_name = basename(entry)
        file_name = file_name.replace('.header', '')
        super().__init__(file_name=file_name, source_names=[entry])

    @property
    def prev(self):
        return '{}.jpg'.format(self._file_id)

    @property
    def thumb(self):
        return '{}_th.jpg'.format(self._file_id)

    def is_valid(self):
        # over-ride self._obs_id dependency
        return True

    def set_destination_uris(self):
        super().set_destination_uris()
        temp = []
        for entry in self._destination_uris:
            temp.append(entry.replace('.header', ''))
        self._destination_uris = temp

    def set_obs_id(self):
        # must use the builder
        self._obs_id = None


class GemProcBuilder(nbc.StorageNameBuilder):
    def __init__(self, metadata_reader):
        super().__init__()
        self._metadata_reader = metadata_reader
        self._logger = getLogger(self.__class__.__name__)

    def build(self, entry):
        self._logger.debug(f'Begin build for {entry}')
        gem_proc_name = GemProcName(entry=entry)
        gem_proc_name.obs_id = self._get_obs_id(gem_proc_name)
        self._logger.debug(f'End build with {gem_proc_name.file_name}')
        return gem_proc_name

    def _get_obs_id(self, gem_proc_name):
        """
        These files are not available from archive.gemini.edu, so
        only ask for their metadata from CADC.
        """
        self._logger.debug(
            f'Begin _get_obs_id for file_name {gem_proc_name.file_name}'
        )
        self._metadata_reader.set(gem_proc_name)
        headers = self._metadata_reader.headers.get(gem_proc_name.file_uri)
        if headers is None:
            raise CadcException(f'No metadata for {gem_proc_name}')
        data_label = get_keyword(headers, 'DATALAB')
        if data_label is not None:
            data_label = repair_data_label(
                gem_proc_name.file_name, data_label
            )
        self._logger.debug(f'End _get_obs_id with {data_label}')
        return data_label
