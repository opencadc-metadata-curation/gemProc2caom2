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

"""
This module implements the ObsBlueprint mapping, as well as the workflow 
entry point that executes the workflow.
"""

import importlib
import logging
import os
import sys
import traceback

from caom2 import Observation, CalibrationLevel, ProductType
from caom2 import Axis, CoordAxis1D, SpectralWCS, CoordFunction1D, RefCoord
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from gem2caom2 import GemName, ARCHIVE, COLLECTION, SCHEME
from gem2caom2 import external_metadata as em


__all__ = ['gem_proc_main_app', 'update', 'APPLICATION', 'GemProcName',
           'to_caom2']


APPLICATION = 'gemProc2caom2'


class GemProcName(GemName):
    """A class to over-ride the gemini schema for the fits files."""

    def __init__(self, file_name):
        super(GemProcName, self).__init__(file_name=file_name)
        self.scheme = 'ad'
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)
    bp.configure_energy_axis(4)
    bp.configure_polarization_axis(5)
    bp.configure_observable_axis(6)

    # processed pipeline, so they're all Derived
    bp.set('CompositeObservation.members', {})
    bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)

    bp.clear('Plane.provenance.lastExecuted')
    bp.add_fits_attribute('Plane.provenance.lastExecuted', 'DATE')
    bp.set('Plane.provenance.name', 'Nifty4Gemini')
    bp.set('Plane.provenance.producer', 'CADC')
    bp.set('Plane.provenance.reference',
           'https://doi.org/10.5281/zenodo.897014')

    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes (an n:n
    relationship between TDM attributes and CAOM attributes). Must have this
    signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = kwargs.get('headers')
    fqn = kwargs.get('fqn')
    uri = kwargs.get('uri')
    gem_proc_name = None
    if uri is not None:
        gem_proc_name = GemProcName(file_name=mc.CaomName(uri).file_name)
    if fqn is not None:
        gem_proc_name = GemProcName(file_name=os.path.basename(fqn))
    if gem_proc_name is None:
        raise mc.CadcException(f'Need one of fqn or uri defined for '
                               f'{observation.observation_id}')

    for plane in observation.planes.values():
        if plane.product_id != gem_proc_name.product_id:
            continue

        cc.update_plane_provenance_list(
            plane, headers, ['FCMB', 'DCMB', 'SHFILE'], COLLECTION,
            _repair_provenance_value, observation.observation_id)

        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                idx = mc.to_int(part.name)
                header = headers[idx]

                for chunk in part.chunks:
                    if chunk.position is None and chunk.naxis is not None:
                        chunk.naxis = None
                    extname = header.get('EXTNAME')
                    if extname == 'SCI':
                        part.product_type = ProductType.SCIENCE
                        _update_flat_energy(
                            chunk, header, observation.observation_id)
                    elif extname == 'DQ':
                        part.product_type = ProductType.AUXILIARY
                    elif extname == 'VAR':
                        part.product_type = ProductType.NOISE

        # update observation members
        # DB 24-06-20
        # At least for the processed flats only the FCOMBINE x FCMB# datasets
        # should be members.  These are the original individual flats that were
        # combined for the processed product.
        provenance_headers = extract_keywords(headers, ['FCMB'])

        def member_filter(plane_uri):
            result = False
            for entry in provenance_headers.values():
                if entry in plane_uri.uri:
                    result = True
            return result
        cc.update_observation_members_filtered(observation, member_filter)

    logging.debug('Done update.')
    return observation


def _update_flat_energy(chunk, header, obs_id):
    logging.debug(f'Begin _update_flat_energy for {obs_id}.')
    axis = Axis(ctype='WAVE', cunit='Angstrom')
    coord_axis_1d = CoordAxis1D(axis)
    ref_coord = RefCoord(pix=header.get('CRPIX1'),
                         val=header.get('CRVAL1'))
    fn = CoordFunction1D(naxis=header.get('NAXIS1'),
                         delta=header.get('CD1_1'),
                         ref_coord=ref_coord)
    coord_axis_1d.function = fn
    energy = SpectralWCS(axis=coord_axis_1d,
                         specsys='TOPOCENT')
    chunk.energy = energy
    chunk.energy_axis = 1
    chunk.naxis = 1
    logging.debug('End _update_flat_energy.')


def extract_keywords(headers, keyword_list):
    """
    :param headers: astropy FITS headers list
    :param keyword_list: list of keyword prefixes to extract.
    :return:
    """
    result = {}
    for entry in keyword_list:
        for header in headers:
            for keyword in header:
                if keyword.startswith(entry):
                    result[keyword] = header.get(keyword)
    return result


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The artifact URIs for the files to be processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(uri):
            accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.local:
        for ii in args.local:
            file_id = mc.StorageName.remove_extensions(os.path.basename(ii))
            file_name = f'{file_id}.fits'
            result.append(GemProcName(file_name=file_name).file_uri)
    elif args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result


def _repair_provenance_value(value, obs_id):
    prov_file_id = GemName.remove_extensions(value)
    prov_obs_id = em.get_gofr().get_obs_id(prov_file_id)
    prov_obs_id = 'None' if prov_obs_id is None else prov_obs_id
    logging.debug(f'End _repair_provenance_value. {prov_obs_id} '
                  f'{prov_file_id}')
    return prov_obs_id, prov_file_id


def to_caom2():
    """This function is called by pipeline execution. It must have this name.
    """
    args = get_gen_proc_arg_parser().parse_args()
    uris = _get_uris(args)
    blueprints = _build_blueprints(uris)
    result = gen_proc(args, blueprints)
    logging.debug(f'Done {APPLICATION} processing.')
    return result
           

def gem_proc_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
