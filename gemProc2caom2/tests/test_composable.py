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

import os
from test_main_app import TEST_DATA_DIR

from mock import Mock, patch

from caom2pipe.astro_composable import make_headers_from_file
from gemProc2caom2 import composable, GemProcName

TEST_OBS_ID = 'GN-2014A-Q-85-16-003-RGN-FLAT'


@patch('gemProc2caom2.composable.ClientCollection')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run(run_mock, access_mock, clients_mock):
    run_mock.return_value = 0
    access_mock.return_value = 'https://localhost:8080'
    test_f_id = 'test_rgn_flat'
    test_f_name = f'{test_f_id}.fits'

    clients_mock.return_value.data_client.get_head.return_value = (
        make_headers_from_file(
            f'{TEST_DATA_DIR}/wrgnN20070626S0205_arc.fits.header'
        )
    )

    getcwd_orig = os.getcwd
    test_fqn = os.path.join(TEST_DATA_DIR, 'run_test')
    os.getcwd = Mock(return_value=test_fqn)
    try:
        # execution
        test_result = composable._run()
        assert test_result == 0, 'successful execution'
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(test_storage, GemProcName), type(test_storage)
        # this is a nonsense observation ID, but it's what the code
        # returns with the file name of test_rgn_flat.fits, and a
        # starting DATALAB value of GN-2006A-SV-131-11-001-ARC
        assert (
            test_storage.obs_id == 'GN-2006A-SV-131-11-001-ARC-RGN-FLAT'
        ), f'wrong obs id {test_storage.obs_id}'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert test_storage.source_names[0] == os.path.join(
            test_fqn, test_f_name
        ), 'wrong source name'
        clients_mock.return_value.data_client.get_head.assert_called_with(
            'cadc:TEST/test_rgn_flat.fits'
        ), 'get_head'
        clients_mock.return_value.data_client.info.assert_called_with(
            'cadc:TEST/test_rgn_flat.fits'
        ), 'info'
        # no repo client calls because they're all masked by the do_one mock
    finally:
        os.getcwd = getcwd_orig
        for entry in [
            'failure_log.txt',
            'rejected.yml',
            'retries.txt',
            'run_test_report.txt',
            'success_log.txt',
        ]:
            entry_fqn = f'{test_fqn}/{entry}'
            if os.path.exists(entry_fqn):
                os.unlink(entry_fqn)
