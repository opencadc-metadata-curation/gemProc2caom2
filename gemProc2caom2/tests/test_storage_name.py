# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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
#  $Revision: 4 $
#
# ***********************************************************************
#

from mock import patch, Mock

from cadcdata import FileInfo
from caom2utils.data_util import (
    get_local_file_headers,
    make_headers_from_string,
)
from caom2pipe.manage_composable import StorageName
from caom2pipe.reader_composable import VaultReader
from gem2caom2.obs_file_relationship import (
    get_prefix,
    get_suffix,
    repair_data_label,
)
from gemProc2caom2 import builder

from test_main_app import TEST_DATA_DIR


def test_builder():

    original_scheme = StorageName.scheme
    original_collection = StorageName.collection

    try:
        StorageName.scheme = builder.CADC_SCHEME
        StorageName.collection = builder.COLLECTION
        test_id = 'ctfbrsnN20140428S0086'
        test_f_name = f'{test_id}.fits'
        mock_client = Mock()
        test_metadata_reader = VaultReader(mock_client)
        from caom2pipe.astro_composable import make_headers_from_file

        test_fqn = f'{TEST_DATA_DIR}/{test_f_name}.header'
        test_headers = make_headers_from_file(test_fqn)
        test_metadata_reader._retrieve_headers = Mock()
        test_metadata_reader._retrieve_headers.return_value = test_headers
        test_metadata_reader._retrieve_file_info = Mock(
            return_value=FileInfo(
                id=test_f_name,
                size=123,
                md5sum='123',
                file_type='application/fits',
                encoding=None,
            )
        )

        test_subject = builder.GemProcBuilder(test_metadata_reader)
        for entry in [
            test_f_name,
            f'vos:goliaths/tests/{test_f_name}',
            f'/tmp/{test_f_name}',
        ]:
            test_sn = test_subject.build(entry)
            assert test_sn.file_uri == f'cadc:GEMINICADC/{test_f_name}'
            assert test_sn.prev == f'{test_id}.jpg'
            assert test_sn.thumb == f'{test_id}_th.jpg'
            assert test_sn.prev_uri == f'cadc:GEMINICADC/{test_id}.jpg'
            assert test_sn.thumb_uri == f'cadc:GEMINICADC/{test_id}_th.jpg'
            assert test_sn.destination_uris == [
                'cadc:GEMINICADC/ctfbrsnN20140428S0086.fits'
            ], f'wrong destination uris for {entry}'
            assert (
                test_sn.obs_id == 'GN-2014A-Q-85-12-002'
            ), f'wrong obs_id for {entry}'
            assert test_sn._source_names[0] == entry, 'wrong source name'
    finally:
        StorageName.scheme = original_scheme
        StorageName.collection = original_collection


def test_repair():
    candidates = {
        'wrgnN20140428S0085_arc': 'GN-2014A-Q-85-12-001-ARC',
        'wrgnN20070116S0165_arc': 'GN-2006B-C-4-29-015-ARC',
    }
    for entry in candidates.keys():
        fqn = f'{TEST_DATA_DIR}/{entry}.fits'
        headers = get_local_file_headers(fqn)
        data_label = headers[0].get('DATALAB')
        assert get_prefix(entry) == 'wrgn', f'{get_prefix(entry)}'
        test_result = get_suffix(entry, data_label)
        assert test_result == ['arc'], f'{test_result}'

        test_result = repair_data_label(entry, data_label)
        assert test_result == candidates.get(
            entry
        ), f'expected {candidates.get(entry)}'
