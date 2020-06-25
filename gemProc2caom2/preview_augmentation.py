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

import logging
import os
import re

import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ZScaleInterval

from caom2 import ProductType, ReleaseType
from caom2pipe import manage_composable as mc
from gem2caom2 import ARCHIVE, GemName


class GemProcPreview(mc.PreviewVisitor):

    def __init__(self, **kwargs):
        super(GemProcPreview, self).__init__(
            ARCHIVE, ReleaseType.DATA, **kwargs)
        self._storage_name = GemName(file_name=self._science_file)
        self._science_fqn = os.path.join(self._working_dir,
                                         self._storage_name.file_name)
        self._preview_fqn = os.path.join(
            self._working_dir,  self._storage_name.prev)
        self._thumb_fqn = os.path.join(
            self._working_dir, self._storage_name.thumb)
        self._logger = logging.getLogger(__name__)

    def generate_plots(self, obs_id):
        self._logger.error(f'Begin generate_plots for {obs_id}')
        count = 0
        hdus = fits.open(self._science_fqn)
        data_label = hdus[0].header.get('DATALAB').upper()
        interval = ZScaleInterval()
        if 'CTFBRSN' in data_label:
            white_light_data = np.flipud(np.median(hdus['SCI'].data, axis=0))
        elif ('FLAT' in data_label or 'ARC' in data_label or
              'RONCHI' in data_label):
            # Stitch together the 29 'SCI' extensions into one array and save.
            hdul = [x for x in hdus if x.name == 'SCI']
            hdul.sort(key=lambda x: int(re.split(r"[\[\]\:\,']+",
                                                 x.header['NSCUTSEC'])[3]))
            temp = np.concatenate([x.data for x in hdul])
            white_light_data = interval(temp)
        elif 'SHIFT' in data_label:
            temp = np.flipud(hdus['SCI'].data)
            white_light_data = interval(temp)
        else:
            return count

        plt.imsave(self._preview_fqn, white_light_data, cmap='inferno')
        count += 1
        self.add_preview(self._storage_name.prev_uri, self._storage_name.prev,
                         ProductType.PREVIEW)
        self.add_to_delete(self._preview_fqn)
        count += self._gen_thumbnail()
        self._logger.error('End generate_plots.')
        return count

    def _gen_thumbnail(self):
        self._logger.debug(f'Generating thumbnail for file '
                           f'{self._science_fqn}.')
        count = 0
        if os.path.exists(self._preview_fqn):
            thumb = image.thumbnail(self._preview_fqn, self._thumb_fqn,
                                    scale=0.25)
            if thumb is not None:
                self.add_preview(self._storage_name.thumb_uri,
                                 self._storage_name.thumb,
                                 ProductType.THUMBNAIL)
                self.add_to_delete(self._thumb_fqn)
                count = 1
        else:
            self._logger.warning(f'Could not find {self._preview_fqn} for '
                                 f'thumbnail generation.')
        return count


def visit(observation, **kwargs):
    previewer = GemProcPreview(**kwargs)
    return previewer.visit(observation, previewer._storage_name)
