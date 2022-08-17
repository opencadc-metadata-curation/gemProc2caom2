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

from os import path

import logging
from astropy.io import fits
from caom2 import (
    Observation,
    DerivedObservation,
    ObservationURI,
    PlaneURI,
    TypedSet,
)
from caom2pipe.client_composable import query_tap_client, repo_get
from caom2pipe.manage_composable import (
    build_uri,
    CadcException,
    CaomName,
    check_param,
    Config,
    get_keyword,
    TaskType,
)
from gem2caom2.gem_name import COLLECTION, GemName, SCHEME
from gem2caom2.obs_file_relationship import repair_data_label


def visit(observation, **kwargs):
    check_param(observation, Observation)
    storage_name = kwargs.get('storage_name')
    if storage_name is None:
        raise CadcException(
            f'Must have a storage_name parameter for provenance_augmentation '
            f'for {observation.observation_id}'
        )
    clients = kwargs.get('clients')
    if clients is None:
        logging.warning(f'Provenance augmentation does not work for SCRAPE.')
        return observation

    working_directory = kwargs.get('working_directory')
    if working_directory is None:
        raise CadcException(f'Must have a working_directory parameter for provenance_augmentation.')

    obs_members = TypedSet(
        ObservationURI,
    )

    for plane in observation.planes.values():
        plane_inputs = TypedSet(
            PlaneURI,
        )
        for artifact in plane.artifacts.values():
            if storage_name.file_uri == artifact.uri:
                _do_provenance(
                    working_directory,
                    storage_name,
                    observation,
                    plane_inputs,
                    obs_members,
                    clients,
                )

        if plane.provenance is not None:
            plane.provenance.inputs.update(plane_inputs)

    if isinstance(observation, DerivedObservation):
        observation.members.update(obs_members)
        if len(observation.members) > 0:
            observable = kwargs.get('observable')
            _do_members_metadata(
                observation,
                clients,
                observation.members,
                observable.metrics,
            )

    logging.info(
        f'Done provenance_augmentation for {observation.observation_id}'
    )
    return observation


def _do_provenance(
    working_directory,
    storage_name,
    observation,
    plane_inputs,
    obs_members,
    clients,
):
    """
    DB 06-08-20
    Looking at the DATALAB values for the test set, these are now set to
    correctly identify the correct observation ID.
    e.g. rnN20140428S0174_dark.fits has DATALAB = GN-2014A-Q-85-16-006 since
    it is NOT derived from multiple observations.  rgnN20140428S0174_dark.fits
    (with the extra ‘g’) has DATALAB = GN-2014A-Q-85-16-006-DARK since it IS a
    new derived observation.   And inputs/members are in the PROVENANCE
    extension.

    DB 07-08-20
    All members + inputs in the extension are plane.provenance.inputs.

    Add the appropriate raw planes of the ‘member’ observations to the list
    of inputs for the derived observations, where 'appropriate' means
    find the ‘appropriate’ filename that identifies the plane of the inputs.
    """
    science_fqn = storage_name.get_file_fqn(working_directory)
    logging.debug(f'Begin _do_provenance for {observation.observation_id} from {science_fqn}')
    count = 0
    hdus = fits.open(science_fqn)
    if 'PROVENANCE' not in hdus:
        logging.warning(f'PROVENANCE extension not found in HDUs for {science_fqn}.')
        return count

    data = hdus['PROVENANCE'].data
    temp = None
    for entry in data.columns:
        if entry.name.startswith('Type'):
            temp = entry.name
            break
    for f_name, f_prov_type in zip(data['Filename'], data[temp]):
        f_id = GemName.remove_extensions(f_name)

        # there are file names in the data
        # some are GEMINICADC
        # some are GEMINI
        # need the observation IDs (DATALAB values) for the file names
        # start with CAOM2 records, and then do a get_head if necessary

        # query caom2 - have access to a tap client
        # write a query where it doesn't matter what the collection is
        # where A.uri LIKE %{f_name}

        qs = f"""
        SELECT O.observationID, A.uri
        FROM caom2.Observation AS O
        JOIN caom2.Plane AS P ON P.obsID = O.obsID
        JOIN caom2.Artifact AS A ON A.planeID = P.planeID
        WHERE A.uri LIKE '%/{f_name}'
        """
        result = query_tap_client(qs, clients.query_client)
        prov_obs_id = None
        # preserve the prov_uri because it identifies the collection of
        # the provenance
        prov_uri = None
        if len(result) > 0:
            for entry in result:
                prov_obs_id = entry['observationID']
                prov_uri = entry['uri']
                logging.info(f'Found observationID {prov_obs_id}, uri {prov_uri} for {f_name}')
        if len(result) > 1:
            raise CadcException(f'Too many rows:\n{result}')

        if prov_obs_id is None:
            # query SI with get_head
            # try both collections
            # start with GEMINICADC because that's quicker, then GEMINI
            for uri_prefix in ['cadc:GEMINICADC', 'gemini:GEMINI']:
                prov_uri = f'{uri_prefix}/{f_name}'
                headers = clients.data_client.get_head(prov_uri)
                if headers is not None:
                    temp = get_keyword(headers, 'DATALAB')
                    if temp is not None:
                        prov_obs_id = repair_data_label(f_name, temp)
                        logging.info(
                            f'Found observationID {prov_obs_id}, uri '
                            f'{prov_uri} for {f_name} in storage.'
                        )
                        break

        if prov_obs_id is not None and prov_uri is not None:
            logging.info(f'Found observation ID {prov_obs_id} for file {f_id} with provenance type {f_prov_type}.')
            collection = prov_uri.split(':')[1].split('/')[0]
            input_obs_uri_str = CaomName.make_obs_uri_from_obs_id(
                collection, prov_obs_id
            )
            input_obs_uri = ObservationURI(input_obs_uri_str)
            plane_uri = PlaneURI.get_plane_uri(input_obs_uri, f_id)
            plane_inputs.add(plane_uri)
            count += 1
            if f_prov_type == 'member' and isinstance(
                observation, DerivedObservation
            ):
                member_obs_uri_str = CaomName.make_obs_uri_from_obs_id(
                    collection, prov_obs_id
                )
                member_obs_uri = ObservationURI(member_obs_uri_str)
                obs_members.add(member_obs_uri)
                count += 1
    hdus.close()
    logging.debug('End _do_provenance.')
    return count


def _do_members_metadata(observation, clients, members, metrics):
    """
    Look up the metadata for the provenance member, and use that Proposal
    and Moving Target information for the co-added products as well.

    DB - 07-08-20

    * Proposal information for the co-added products (the ones with ‘g’).
    Perhaps the DATALAB isn’t being parsed properly for the original proposal
    ID?

    * ‘Moving target’ could in principal be determined by looking for
    NON_SIDEREAL in the ‘types’ returned by jsonsummary for the unprocessed
    dataset as in gem2caom2.  e.g.
    https://archive.gemini.edu/jsonsummary/GN-2014A-Q-85-12-002.  e.g The
    Titan observation should be shown as being a moving target whereas the
    other object exposure is not.
    """
    logging.debug(
        f'Begin _do_members_metadata for ' f'{observation.observation_id}.'
    )
    prov_obs_uri = None
    for entry in members:
        # use the first one for querying
        prov_obs_uri = entry
        break
    if clients.metadata_client is not None:
        prov_obs = repo_get(
            clients.metadata_client,
            prov_obs_uri.collection,
            prov_obs_uri.observation_id,
            metrics,
        )
        if prov_obs is not None:
            observation.proposal = prov_obs.proposal
            if prov_obs.target is not None:
                observation.target.moving = prov_obs.target.moving

    logging.debug(f'Done _do_members_metadata.')
