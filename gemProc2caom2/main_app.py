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

DB 21-07-20
CTFBRSN files:

The first extension is the processed NIFS data cube. Axes 1 and 2 are
spatial axes.  Axis 3 is spectral.

The second extension is a binary table.  No WCS information.

3rd extension is a bad pixel map. WCS info should be ignored for this.

processed pipeline, so even though some might think they're all Derived:
DB/NC 06-07-20

How to know if a NIFS arc file is a unique observation:
If it’s a co-add of more than one individual observation then it is a
new composite (maybe ‘derived’ is the new term) observation. If instead
the processed arc is produced from just a single unprocessed observation
then it is NOT a new observation but just a processed version of the
original unprocessed dataset, or a new plane.

The G in WRGN means the file passed through the Co-adding task. It's
just that sometimes only one file is passed to and it's not really
co-adding anything. It would make more sense to only add the G prefix
if there was actually Co-adding happening.

Modify the Nifty code to only add that G for co-added exposures. Then
the datalabel for a co-added processed arc would look like
GN-2014A-Q-85-12-001-WRGN-ARC and a non-co-added processed arc would
look like GN-2014A-Q-85-12-001-WRN-ARC

DB 06-08-20
If only one ‘member’ is present then it’s a new CAL=2 plane for the
    original observation given by the member name.  If there is > 1 ‘member’
then it’s a new derived observation.

"""

import importlib
import logging
import os
import sys
import traceback

from cadctap import CadcTapClient
from cadcutils import net
from caom2 import Observation, CalibrationLevel, ProductType, TemporalWCS
from caom2 import Axis, CoordAxis1D, SpectralWCS, CoordFunction1D, RefCoord
from caom2 import CoordError, ObservationIntentType, SimpleObservation
from caom2 import Algorithm, DataProductType, DerivedObservation, Provenance
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import WcsParser
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from gem2caom2 import ARCHIVE, COLLECTION, external_metadata
from gem2caom2 import gem_name, obs_file_relationship


__all__ = ['gem_proc_main_app', 'update', 'APPLICATION', 'GemProcName',
           'to_caom2']


APPLICATION = 'gemProc2caom2'


class GemProcName(mc.StorageName):

    def __init__(self, file_name, entry):
        super(GemProcName, self).__init__(fname_on_disk=file_name,
                                          archive=ARCHIVE,
                                          compression='',
                                          entry=entry)
        self._file_name = file_name
        self._file_id = gem_name.GemName.remove_extensions(file_name)
        self._obs_id = self.get_obs_id()
        self.scheme = 'cadc'
        self.archive = 'GEMINI'
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    def get_obs_id(self):
        obs_id = external_metadata.get_obs_id_from_headers(self._file_id)
        if obs_id is None:
            config = mc.Config()
            config.get_executors()
            subject = mc.define_subject(config)
            tap_client = CadcTapClient(subject=subject,
                                       resource_id=config.tap_id)
            obs_id = external_metadata.get_obs_id_from_cadc(
                self._file_id, tap_client)
        return obs_id

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
        return self._file_id

    @property
    def thumb(self):
        return '{}_th.jpg'.format(self._file_id)


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1, 2))
    bp.configure_energy_axis(3)
    bp.configure_polarization_axis(5)
    bp.configure_observable_axis(6)

    # assume all observations are Derived, set to Simple if required
    bp.set('DerivedObservation.members', {})

    bp.set('Observation.observationID', '_get_obs_id(header)')
    bp.set('Observation.intent', '_get_obs_intent(uri)')

    # TODO - this is not correct, as it will over-write the value in the
    # existing gem2caom2 version
    bp.clear('Observation.algorithm.name')
    bp.add_fits_attribute('Observation.algorithm.name', 'SOFTWARE')

    # DB 12-04-21
    # Proposal information comes from GEMPRGID and the details should be
    # retrieved from archive.gemini.edu
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'GEMPRGID')

    # DB 07-08-20
    # target.type for all NIFS science products, at least, could be set to
    # ‘object’.   (i.e. NIFS is never used for wider ‘field’ observations)
    bp.set('Observation.target.type', 'object')
    bp.set('Observation.telescope.geoLocationX', '_get_telescope_x(uri)')
    bp.set('Observation.telescope.geoLocationY', '_get_telescope_y(uri)')
    bp.set('Observation.telescope.geoLocationZ', '_get_telescope_z(uri)')

    bp.set('Plane.calibrationLevel', '_get_plane_calibration_level(header)')
    bp.set('Plane.dataProductType', '_get_plane_data_product_type(header)')
    bp.clear('Plane.provenance.lastExecuted')
    bp.add_fits_attribute('Plane.provenance.lastExecuted', 'DATE')
    bp.clear('Plane.provenance.name')
    bp.add_fits_attribute('Plane.provenance.name', 'SOFTWARE')
    bp.set('Plane.provenance.producer', 'CADC')
    bp.clear('Plane.provenance.reference')
    bp.add_fits_attribute('Plane.provenance.reference', 'SOFT_DOI')
    bp.clear('Plane.provenance.version')
    bp.add_fits_attribute('Plane.provenance.version', 'SOFT_VER')

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
        temp = mc.CaomName(uri).file_name
        gem_proc_name = GemProcName(file_name=temp, entry=temp)
    if fqn is not None:
        temp = os.path.basename(fqn)
        gem_proc_name = GemProcName(file_name=temp, entry=temp)
    if gem_proc_name is None:
        raise mc.CadcException(f'Need one of fqn or uri defined for '
                               f'{observation.observation_id}')

    for plane in observation.planes.values():
        if plane.product_id != gem_proc_name.product_id:
            continue

        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                idx = mc.to_int(part.name)
                header = headers[idx]
                extname = header.get('EXTNAME')
                # DB 22-07-20
                # There are a few other EXTNAME values to look at for
                # part.ProductType.   MDF values would be ‘AUXILIARY’.  The
                # ones currently called “CAL” are likely best set to ‘INFO’
                # since it contains info about datasets used to produce the
                # product.
                #
                # DB 07-08-20
                # EXTNAME  in (‘DQ’, ‘VAR’) should both have
                # ProductType.NOISE.   ‘CAL’ should no longer exist - it’s now
                # BPM. Default type is 'AUXILIARY', 'SCI' is type 'SCIENCE'
                if extname == 'SCI':
                    part.product_type = ProductType.SCIENCE
                elif extname in ['DQ', 'VAR']:
                    part.product_type = ProductType.NOISE
                else:
                    part.product_type = ProductType.AUXILIARY

                if (part.product_type in
                        [ProductType.SCIENCE, ProductType.INFO]):
                    for chunk in part.chunks:
                        filter_name = headers[0].get('FILTER').split('_')[0]
                        _update_energy(
                            chunk, headers[idx], filter_name,
                            observation.observation_id)
                        _update_time(part, chunk, headers[0],
                                     observation.observation_id)
                        if part.product_type == ProductType.SCIENCE:
                            _update_spatial_wcs(part, chunk, headers,
                                                observation.observation_id)
                            chunk.naxis = header.get('NAXIS')
                            if (chunk.position is None and
                                    chunk.naxis is not None):
                                chunk.naxis = None

                        if (chunk.time is not None and
                                chunk.time.axis is not None and
                                chunk.time.axis.function is not None and
                                chunk.time.axis.function.delta == 1.0):
                            # these are the default values, and they make
                            # the time range start in 1858
                            chunk.time = None
                else:
                    # DB 21-07-20
                    # ignore WCS information unless product type == SCIENCE
                    while len(part.chunks) > 0:
                        del part.chunks[-1]

        if isinstance(observation, DerivedObservation):
            if plane.provenance is None:
                plane.provenance = Provenance(name='TBD')
            subject = net.Subject(certificate='/usr/src/app/cadcproxy.pem')
            tap_client = CadcTapClient(
                    subject, resource_id='ivo://cadc.nrc.ca/ams/gemini')

            def _repair_provenance_value(imcmb_value, obs_id):
                logging.debug(f'Being _repair_provenance_value for {obs_id}.')
                prov_file_id = gem_name.GemName.remove_extensions(imcmb_value)
                prov_obs_id = external_metadata.get_obs_id_from_cadc(
                    prov_file_id, tap_client)
                logging.debug(f'End _repair_provenance_value {prov_obs_id} '
                              f'{prov_file_id}')
                return prov_obs_id, prov_file_id

            cc.update_plane_provenance_list(
                    plane, headers, 
                    ['IMCMB', 'SKY', 'FLATIM', 'DARKIM', 'BPMIMG'], 
                    COLLECTION, _repair_provenance_value, 
                    observation.observation_id)
            cc.build_temporal_wcs_bounds(tap_client, headers[0],
                    ['IMCMB', 'SKY', 'FLATIM', 'DARKIM', 'BPMIMG'],
                    COLLECTION)

    if (observation.proposal is not None and
            observation.proposal.id is not None and
            observation.proposal.pi_name is None):
        program = external_metadata.get_pi_metadata(observation.proposal.id)
        if program is not None:
            observation.proposal.pi_name = program.get('pi_name')
            observation.proposal.title = program.get('title')

    if isinstance(observation, SimpleObservation):
        # undo the observation-level metadata modifications for updated
        # Gemini records
        observation.algorithm = Algorithm(name='exposure')
    else:
        cc.update_observation_members(observation)
    logging.debug('Done update.')
    return observation


def _get_obs_id(header):
    return header.get('DATALAB')


def _get_obs_intent(uri):
    # DB 07-08-20
    # The co-added products (the ones with ‘g’ in the file name?) the Intent
    # should be set to calibration.
    ignore_scheme, ignore_path, file_name = mc.decompose_uri(uri)
    prefix = obs_file_relationship.get_prefix(file_name)
    result = ObservationIntentType.SCIENCE
    if 'g' in prefix:
        result = ObservationIntentType.CALIBRATION
    return result


def _get_plane_calibration_level(header):
    result = CalibrationLevel.RAW_STANDARD
    for keyword in ['IMCMB', 'SKY', 'FLATIM', 'DARKIM', 'BPMIMG']:
        for key in header:
            if key.startswith(keyword):
                result = CalibrationLevel.CALIBRATED
                break
    return result


def _get_plane_data_product_type(header):
    # DB 19-08-20
    # Any with OBSTYPE of OBJECT should be set to ‘cube’ now. Products with
    # OBSTYPE’s of ARC/RONCHI/FLAT/DARK (i.e. all non-OBJECT) should be set to
    # ‘spectrum’,  since one axis is the wavelengths axis, and more consistent
    # with the final use of the data (and this is what data product type is set
    # to for the raw data).
    result = DataProductType.SPECTRUM
    obs_type = header.get('OBSTYPE')
    if obs_type == 'OBJECT':
        instrument = header.get('INSTRUME')
        if instrument == 'NIRI':
            result = DataProductType.IMAGE
        else:
            result = DataProductType.CUBE
    return result


def _get_telescope_x(uri):
    x, ignore_y, ignore_z = _get_telescope()
    return x


def _get_telescope_y(uri):
    ignore_x, y, ignore_z = _get_telescope()
    return y


def _get_telescope_z(uri):
    ignore_x, ignore_y, z = _get_telescope()
    return z


def _get_telescope():
    # DB 19-08-20
    # NIFS only ever on Gemini North
    x, y, z = ac.get_location(19.823806, -155.46906, 4213.0)
    return x, y, z


def _update_energy(chunk, header, filter_name, obs_id):
    logging.debug(f'Begin _update_energy for {obs_id}.')
    # because the type for the axes are 'LINEAR', which isn't an energy type,
    # so can't use the WcsParser from caom2utils.
    disp_axis = header.get('DISPAXIS')
    naxis = header.get('NAXIS')
    if disp_axis is not None and naxis is not None and disp_axis <= naxis:
        axis = Axis(ctype='WAVE', cunit='Angstrom')
        coord_axis_1d = CoordAxis1D(axis)
        ref_coord = RefCoord(pix=header.get(f'CRPIX{disp_axis}'),
                             val=header.get(f'CRVAL{disp_axis}'))
        fn = CoordFunction1D(naxis=header.get(f'NAXIS{disp_axis}'),
                             delta=header.get(f'CD{disp_axis}_{disp_axis}'),
                             ref_coord=ref_coord)
        coord_axis_1d.function = fn
        energy = SpectralWCS(axis=coord_axis_1d,
                             specsys='TOPOCENT')
        energy.bandpass_name = filter_name
        # DB 07-08-20
        # I think the best we can do is assume that a resolution element is 2
        # pixels wide. So resolving power is the absolute value of
        # approximately CRVAL3/(2 * CD3_3)
        energy.resolving_power = abs(header.get(f'CRVAL{disp_axis}') / (
                    2 * header.get(f'CD{disp_axis}_{disp_axis}')))
        chunk.energy = energy
        chunk.energy_axis = disp_axis
    logging.debug('End _update_energy.')


def _update_spatial_wcs(part, chunk, headers, obs_id):
    logging.debug(f'Begin _update_spatial_wcs for {obs_id}')
    if part.name == '1':
        # DB/NC 22-07-20
        #
        # Use the NAXIS* values from the 'SCI' extension.
        #
        # CD*:  That should give approximately correct pixel scales along the
        # axes in degrees/pixel.
        #
        # The hard coded values for CRPIX assume the RA/DEC values are at the
        # centre of the 60/62 array of pixels in the final product.  The
        # primary header CRPIX values would have referred to the original
        # ~2000 x ~2000 pixel detector array I’m guessing.

        header = headers[0]
        idx = mc.to_int(part.name)
        header['NAXIS1'] = headers[idx].get('NAXIS1')
        header['NAXIS2'] = headers[idx].get('NAXIS2')
        header['CD1_1'] = 3.0 / (header.get('NAXIS1') * 3600.0)
        header['CD2_2'] = 3.0 / (header.get('NAXIS2') * 3600.0)
        header['CD1_2'] = 0.0
        header['CD2_1'] = 0.0
        header['CRPIX1'] = header.get('NAXIS1') / 2.0
        header['CRPIX2'] = header.get('NAXIS2') / 2.0
        wcs_parser = WcsParser(header, obs_id, 0)
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
    logging.debug(f'End _update_spatial_wcs')


def _update_time(part, chunk, header, obs_id):
    logging.debug(f'Begin _update_time for {obs_id} part {part.name}.')

    # DB 02-07-20
    # time metadata comes from MJD_OBS and EXPTIME, it's not
    # an axis requiring cutout support
    exp_time = header.get('EXPTIME')
    mjd_obs = header.get('MJD_OBS')
    if exp_time is None or mjd_obs is None:
        chunk.time = None
    else:
        if chunk.time is None:
            coord_error = CoordError(syser=1e-07, rnder=1e-07)
            time_axis = CoordAxis1D(axis=Axis('TIME', 'd'), error=coord_error)
            chunk.time = TemporalWCS(axis=time_axis,
                                     timesys='UTC')
        ref_coord = RefCoord(pix=0.5, val=mjd_obs)
        chunk.time.axis.function = CoordFunction1D(
            naxis=1,
            delta=mc.convert_to_days(exp_time),
            ref_coord=ref_coord)
        chunk.time.exposure = exp_time
        chunk.time.resolution = mc.convert_to_days(exp_time)
    logging.debug(f'End _update_time.')


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
        blueprint = ObsBlueprint(module=module, update=True)
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
            result.append(GemProcName(file_name=file_name,
                                      entry=file_name).file_uri)
    elif args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result


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
