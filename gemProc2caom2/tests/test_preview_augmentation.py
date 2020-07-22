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

from mock import patch

from caom2pipe import manage_composable as mc
from gemProc2caom2 import preview_augmentation
import test_main_app

REJECTED_FILE = os.path.join(test_main_app.TEST_DATA_DIR, 'rejected.yml')
TEST_FILES_DIR = '/test_files'


@patch('caom2pipe.manage_composable.query_tap_client')
@patch('caom2utils.fits2caom2.CadcDataClient')
def test_preview_augmentation(data_client_mock, tap_mock):
    tap_mock.side_effect = test_main_app._tap_mock
    data_client_mock.return_value.get_file_info.side_effect = \
        test_main_app._get_file_info

    test_f_id = 'rgnN20140428S0171_flat'
    test_f_name = f'{test_f_id}.fits'
    test_obs = mc.read_obs_from_file(
        f'{test_main_app.TEST_DATA_DIR}/'
        f'{test_f_id}.expected.xml')

    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    kwargs = {'working_directory': TEST_FILES_DIR,
              'cadc_client': None,
              'stream': 'stream',
              'observable': test_observable,
              'science_file': test_f_name}

    try:
        from datetime import datetime
        import logging
        start_ts = datetime.utcnow().timestamp()
        test_result = preview_augmentation.visit(test_obs, **kwargs)
        end_ts = datetime.utcnow().timestamp()
        logging.error(f'{test_f_name} execution time {end_ts - start_ts}')
    except Exception as e:
        import logging
        logging.error(e)
        import traceback
        logging.error(traceback.format_exc())
        assert False

    assert test_result is not None, 'expect a result'
    assert test_result.get('artifacts') == 2, 'wrong result'
