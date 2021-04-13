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
#  : 4 $
#
# ***********************************************************************
#

import os
from mock import Mock, patch

from caom2pipe import manage_composable as mc
from gemProc2caom2 import provenance_augmentation
import test_main_app

REJECTED_FILE = os.path.join(test_main_app.TEST_DATA_DIR, 'rejected.yml')


def pytest_generate_tests(metafunc):
    fqn_list = [f'{test_main_app.TEST_DATA_DIR}/'
                f'rnN20140428S0181_ronchi.expected.xml']
    metafunc.parametrize('test_fqn', fqn_list)


@patch('caom2pipe.manage_composable.repo_get')
@patch('gem2caom2.external_metadata.get_obs_id_from_cadc')
def test_provenance_augmentation(obs_id_mock, repo_get_mock, test_fqn):
    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    repo_get_mock.side_effect = _repo_get_mock
    obs_id_mock.side_effect = _get_obs_id_mock
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_main_app.TEST_DATA_DIR)
    try:
        test_obs = mc.read_obs_from_file(test_fqn)
        assert not test_obs.target.moving, 'initial conditions moving target'
        kwargs = {'science_file': os.path.basename(test_fqn).replace(
                                    '.expected.xml', '.fits'),
                  'working_directory': test_main_app.TEST_DATA_DIR,
                  'observable': test_observable,
                  'caom_repo_client': Mock()}
        test_result = provenance_augmentation.visit(test_obs, **kwargs)
        assert test_result is not None, 'expect a result'
        assert test_result.get('provenance') == 2, 'wrong result'
        assert len(test_obs.members) == 1, 'wrong membership'
        assert test_obs.target.moving, 'should be changed'
    finally:
        os.getcwd = getcwd_orig


def _get_obs_id_mock(f_id, ignore):
    lookup = {'N20140428S0181': 'GN-2014A-Q-85-16-013'}
    return lookup.get(f_id)


def _repo_get_mock(ignore1, ignore2, ignore3, ignore4):
    return mc.read_obs_from_file(
        f'{test_main_app.TEST_DATA_DIR}/GN-2014A-Q-85-16-013.xml')
