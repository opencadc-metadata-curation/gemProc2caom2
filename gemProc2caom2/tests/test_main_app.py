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

from mock import patch, Mock

from astropy.table import Table
from cadcdata import FileInfo
from gem2caom2 import external_metadata as em
from gemProc2caom2 import main_app, GemProcName, APPLICATION, COLLECTION
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

import glob
import logging
import os
import sys
import traceback

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


@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('gem2caom2.program_metadata.get_pi_metadata')
@patch('caom2pipe.client_composable.query_tap_client')
@patch('caom2utils.data_util.StorageClientWrapper')
@patch('gem2caom2.external_metadata.defining_metadata_finder')
def test_main_app(
    obs_id_mock,
    data_client_mock,
    tap_mock,
    get_pi_mock,
    local_headers_mock,
    test_name,
):
    obs_id_mock.return_value.get.side_effect = _obs_id_mock
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DATA_DIR)
    tap_mock.side_effect = _tap_mock
    get_pi_mock.side_effect = _get_pi_mock
    # during operation, want to use astropy on FITS files
    # but during testing want to use headers and built-in Python file
    # operations
    local_headers_mock.side_effect = ac.make_headers_from_file
    try:
        test_config = mc.Config()
        test_config.task_types = [mc.TaskType.SCRAPE]
        test_config.use_local_files = True
        basename = os.path.basename(test_name)
        file_name = basename.replace('.header', '')
        gem_name = GemProcName(entry=file_name)
        obs_path = f'{TEST_DATA_DIR}/{gem_name.file_id}.expected.xml'
        input_file = f'{TEST_DATA_DIR}/{gem_name.file_id}.in.xml'
        output_file = f'{TEST_DATA_DIR}/{basename}.actual.xml'

        if os.path.exists(output_file):
            os.unlink(output_file)

        local = _get_local(basename)

        data_client_mock.return_value.info.side_effect = _get_file_info

        sys.argv = (
            f'{APPLICATION} --no_validate '
            f'--local {local} -i {input_file} -o '
            f'{output_file} --plugin {PLUGIN} --module {PLUGIN} --lineage '
            f'{_get_lineage(gem_name)}'
        ).split()
        print(sys.argv)
        try:
            main_app.to_caom2()
        except Exception as e:
            logging.error(traceback.format_exc())
            raise e

        compare_result = mc.compare_observations(output_file, obs_path)
        if compare_result is not None:
            raise AssertionError(compare_result)
        # assert False  # cause I want to see logging messages
    finally:
        os.getcwd = getcwd_orig


def _get_file_info(uri):
    return FileInfo(id=uri, file_type='application/fits')


def _get_lineage(blank_name):
    result = mc.get_lineage(
        COLLECTION, blank_name.product_id, f'{blank_name.file_name}', 'cadc'
    )
    return result


def _get_local(obs_id):
    return f'{TEST_DATA_DIR}/{obs_id}'


def _tap_mock(query_string, mock_tap_client):
    if 'observationID' in query_string:
        return Table.read(
            f'observationID,lastModified\n'
            f'GN-2014A-Q-85-16-003-RGN-FLAT,'
            f'2020-02-25T20:36:31.230\n'.split('\n'),
            format='csv',
        )
    else:
        return Table.read(
            'val,delta,cunit,naxis\n'
            '57389.66314699074,0.000115798611111111,d,1\n'.split('\n'),
            format='csv',
        )


def _obs_id_mock(uri):
    ign1, ign2, f_name = mc.decompose_uri(uri)
    result = None
    if f_name in OID_LOOKUP:
        result = em.DefiningMetadata('GNIRS', OID_LOOKUP.get(f_name))
    return result


def _get_pi_mock(ignore):
    return {
        'pi_name': 'Principle Investigator',
        'title': 'The Title For A NIFS Science Program',
    }
