# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from tempfile import TemporaryDirectory
from mock import patch, Mock

from caom2pipe.manage_composable import StorageName
from gemProc2caom2 import fits2caom2_augmentation, GemProcBuilder
from caom2.diff import get_differences
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from gemProc2caom2.builder import CADC_SCHEME, COLLECTION

import glob
import logging
import os

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')

OID_LOOKUP = {
    'N20140428S0179': 'GN-2014A-Q-85-16-011',
    'N20140428S0180': 'GN-2014A-Q-85-16-012',
    'N20140428S0085': 'GN-2014A-Q-85-12-001',
    'N20130530S0250': 'GN-2013A-Q-62-54-001',
    'N20130530S0258': 'GN-2013A-Q-62-55-005',
    'N20130530S0367': 'GN-2013A-Q-62-59-011',
    'N20130530S0368': 'GN-2013A-Q-62-59-012',
    'N20130530S0265': 'GN-2013A-Q-62-62-002',
    'N20191027S0342': 'GN-2019B-FT-101-34-012',
    'N20191027S0341': 'GN-2019B-FT-101-34-011',
    'N20191027S0083': 'GN-2019B-FT-101-31-001',
    'rgnN20140428S0177': 'GN-2014A-Q-85-16-009',
    'rgnN20140428S0171': 'GN-2014A-Q-85-16-003',
    'rgnN20191027S0336': 'GN-2019B-FT-101-34-006',
    'rgnN20191027S0331': 'GN-2019B-FT-101-34-001',
    'N20130530S0363': 'GN-2013A-Q-62-59-007',
    'N20130530S0364': 'GN-2013A-Q-62-59-008',
    'N20130530S0365': 'GN-2013A-Q-62-59-009',
    'N20130530S0366': 'GN-2013A-Q-62-59-010',
    'N20130530S0362': 'GN-2013A-Q-62-59-006',
    'N20130530S0369': 'GN-2013A-Q-62-59-013',
    'N20130530S0370': 'GN-2013A-Q-62-59-014',
    'N20140428S0174': 'GN-2014A-Q-85-16-006',
    'N20140428S0177': 'GN-2014A-Q-85-16-009',
    'N20140428S0181': 'GN-2014A-Q-85-16-013',
    'N20140428S0178': 'GN-2014A-Q-85-16-010',
    'N20140428S0182': 'GN-2014A-Q-85-16-014',
    'N20140428S0175': 'GN-2014A-Q-85-16-007',
    'N20140428S0176': 'GN-2014A-Q-85-16-008',
    'rnN20191231S0495_flat': 'GN-2019B-Q-303-138-004',
    'rnN20191231S0501_dark': 'GN-2019B-Q-303-138-010',
    'N20191231S0505': 'GN-2019B-Q-303-138-014',
}


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{TEST_DATA_DIR}/*.fits.header')
    metafunc.parametrize('test_name', obs_id_list)


@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
def test_visitor(
    access_url,
    get_pi_mock,
    test_name,
):

    access_url.return_value = 'https://localhost:8080'
    get_pi_mock.side_effect = _get_pi_mock

    original_collection = StorageName.collection
    original_scheme = StorageName.scheme
    try:
        StorageName.collection = COLLECTION
        StorageName.scheme = CADC_SCHEME
        with TemporaryDirectory() as tmp_dir_name:
            test_config = mc.Config()
            test_config.task_types = [mc.TaskType.SCRAPE]
            test_config.use_local_files = True
            test_config.data_sources = [TEST_DATA_DIR]
            test_config.working_directory = tmp_dir_name
            test_config.proxy_fqn = f'{tmp_dir_name}/test_proxy.pem'

            with open(test_config.proxy_fqn, 'w') as f:
                f.write('test content')

            metadata_reader = rdc.FileMetadataReader()
            metadata_reader._retrieve_headers = ac.make_headers_from_file
            builder = GemProcBuilder(metadata_reader)
            storage_name = builder.build(test_name)

            kwargs = {
                'storage_name': storage_name,
                'metadata_reader': metadata_reader,
            }
            # logging.getLogger('').setLevel(logging.DEBUG)
            logging.getLogger('FileMetadataReader').setLevel(logging.DEBUG)
            observation = None
            observation = fits2caom2_augmentation.visit(observation, **kwargs)

            expected_fqn = (
                f'{TEST_DATA_DIR}/{storage_name.file_id}.expected.xml'
            )
            actual_fqn = expected_fqn.replace('expected', 'actual')
            mc.write_obs_to_file(observation, actual_fqn)
            expected = mc.read_obs_from_file(expected_fqn)
            compare_result = get_differences(expected, observation)
            if compare_result is not None:
                actual_fqn = expected_fqn.replace('expected', 'actual')
                mc.write_obs_to_file(observation, actual_fqn)
                compare_text = '\n'.join([r for r in compare_result])
                msg = (
                    f'Differences found in observation {expected.observation_id}\n'
                    f'{compare_text}'
                )
                raise AssertionError(msg)
    finally:
        StorageName.scheme = original_scheme
        StorageName.collection = original_collection


def _get_pi_mock(ignore):
    return {
        'pi_name': 'Principle Investigator',
        'title': 'The Title For A NIFS Science Program',
    }
