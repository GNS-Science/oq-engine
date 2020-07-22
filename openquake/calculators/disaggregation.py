# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2015-2020 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Disaggregation calculator core functionality
"""
import logging
import operator
import numpy

from openquake.baselib import parallel, hdf5
from openquake.baselib.general import (
    AccumDict, get_array_nbytes, humansize, pprod, agg_probs, block_splitter)
from openquake.baselib.python3compat import encode
from openquake.hazardlib import stats
from openquake.hazardlib.calc import disagg
from openquake.hazardlib.imt import from_string
from openquake.hazardlib.gsim.base import ContextMaker
from openquake.hazardlib.contexts import RuptureContext
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.source.rupture import RuptureProxy
from openquake.commonlib import util
from openquake.calculators import getters
from openquake.calculators import base

POE_TOO_BIG = '''\
Site #%d: you are trying to disaggregate for poe=%s.
However the source model produces at most probabilities
of %.7f for rlz=#%d, IMT=%s.
The disaggregation PoE is too big or your model is wrong,
producing too small PoEs.'''
U8 = numpy.uint8
U16 = numpy.uint16
U32 = numpy.uint32
F32 = numpy.float32


def _check_curves(sid, rlzs, curves, imtls, poes_disagg):
    # there may be sites where the sources are too small to produce
    # an effect at the given poes_disagg
    for rlz, curve in zip(rlzs, curves):
        for imt in imtls:
            max_poe = curve[imt].max()
            for poe in poes_disagg:
                if poe > max_poe:
                    logging.warning(POE_TOO_BIG, sid, poe, max_poe, rlz, imt)
                    return True


def _matrix(matrices, num_trts, num_mag_bins):
    # convert a dict trti, magi -> matrix into a single matrix
    trti, magi = next(iter(matrices))
    mat = numpy.zeros((num_trts, num_mag_bins) + matrices[trti, magi].shape)
    for trti, magi in matrices:
        mat[trti, magi] = matrices[trti, magi]
    return mat


def _iml4(rlzs, iml_disagg, imtls, poes_disagg, curves):
    # an ArrayWrapper of shape (N, M, P, Z)
    N, Z = rlzs.shape
    P = len(poes_disagg)
    M = len(imtls)
    arr = numpy.empty((N, M, P, Z))
    for m, imt in enumerate(imtls):
        for (s, z), rlz in numpy.ndenumerate(rlzs):
            curve = curves[s][z]
            if poes_disagg == (None,):
                arr[s, m, 0, z] = imtls[imt]
            elif curve:
                poes = curve[imt][::-1]
                imls = imtls[imt][::-1]
                arr[s, m, :, z] = numpy.interp(poes_disagg, poes, imls)
    return hdf5.ArrayWrapper(arr, {'rlzs': rlzs})


def _prepare_ctxs(dstore, rctx, magstr, cmaker):
    sitecol = dstore['sitecol']
    grp = {n: d[rctx['idx']] for n, d in dstore['mag_%s' % magstr].items()
           if n.endswith('_')}
    ctxs = []
    for u, rec in enumerate(rctx):
        ctx = RuptureContext()
        for par in rctx.dtype.names:
            setattr(ctx, par, rec[par])
        for par in grp:
            setattr(ctx, par[:-1], grp[par][u])
        for par in cmaker.REQUIRES_SITES_PARAMETERS:
            setattr(ctx, par, sitecol[par][ctx.sids])
        ctx.idx = {sid: idx for idx, sid in enumerate(ctx.sids)}
        ctxs.append(ctx)
    # sorting for debugging convenience
    ctxs.sort(key=lambda ctx: ctx.occurrence_rate)
    close_ctxs = [[] for sid in sitecol.sids]
    for ctx in ctxs:
        for sid in ctx.idx:
            close_ctxs[sid].append(ctx)
    return ctxs, close_ctxs


def output(mat6):
    """
    :param mat6: a 6D matrix with axis (D, Lo, La, E, P, Z)
    :returns: two matrices of shape (D, E, P, Z) and (Lo, La, P, Z)
    """
    return pprod(mat6, axis=(1, 2)), pprod(mat6, axis=(0, 3))


def compute_disagg(dstore, rctx, cmaker, iml4, trti, magstr, bin_edges, oq,
                   monitor):
    # see https://bugs.launchpad.net/oq-engine/+bug/1279247 for an explanation
    # of the algorithm used
    """
    :param dstore:
        a DataStore instance
    :param rctx:
        an array of rupture parameters
    :param cmaker:
        a :class:`openquake.hazardlib.gsim.base.ContextMaker` instance
    :param iml4:
        an ArrayWrapper of shape (N, M, P, Z)
    :param trti:
        tectonic region type index
    :param magi:
        magnitude bin index
    :param bin_egdes:
        a quartet (dist_edges, lon_edges, lat_edges, eps_edges)
    :param monitor:
        monitor of the currently running job
    :returns:
        a dictionary sid, imti -> 6D-array
    """
    RuptureContext.temporal_occurrence_model = PoissonTOM(
        oq.investigation_time)
    with monitor('reading rupdata', measuremem=True):
        dstore.open('r')
        ctxs, close_ctxs = _prepare_ctxs(dstore, rctx, magstr, cmaker)  # fast

    magi = numpy.searchsorted(bin_edges[0], float(magstr)) - 1
    if magi == -1:  # when the magnitude is on the edge
        magi = 0
    dis_mon = monitor('disaggregate', measuremem=False)
    ms_mon = monitor('disagg mean_std', measuremem=True)
    N, M, P, Z = iml4.shape
    g_by_z = numpy.zeros((N, Z), numpy.uint8)
    for g, rlzs in enumerate(cmaker.gsims.values()):
        for (s, z), r in numpy.ndenumerate(iml4.rlzs):
            if r in rlzs:
                g_by_z[s, z] = g
    eps3 = disagg._eps3(cmaker.trunclevel, oq.num_epsilon_bins)
    res = {'trti': trti, 'magi': magi}
    imts = [from_string(im) for im in oq.imtls]
    with ms_mon:
        # compute mean and std for a single IMT to save memory
        # the size is N * U * G * 8 bytes
        disagg.set_mean_std(ctxs, imts, cmaker.gsims)

    # disaggregate by site, IMT
    for s, iml3 in enumerate(iml4):
        if not close_ctxs[s]:
            continue
        # dist_bins, lon_bins, lat_bins, eps_bins
        bins = (bin_edges[1], bin_edges[2][s], bin_edges[3][s],
                bin_edges[4])
        iml2 = dict(zip(imts, iml3))
        with dis_mon:
            # 7D-matrix #distbins, #lonbins, #latbins, #epsbins, M, P, Z
            matrix = disagg.disaggregate(
                close_ctxs[s], g_by_z[s], iml2, eps3, s, bins)  # 7D-matrix
            for m in range(M):
                mat6 = matrix[..., m, :, :]
                if mat6.any():
                    res[s, m] = output(mat6)
    return res
    # NB: compressing the results is not worth it since the aggregation of
    # the matrices is fast and the data are not queuing up


def get_outputs_size(shapedic, disagg_outputs):
    """
    :returns: the total size of the outputs
    """
    tot = AccumDict(accum=0)
    for out in disagg_outputs:
        tot[out] = 8
        for key in out.lower().split('_'):
            tot[out] *= shapedic[key]
    return tot * shapedic['N'] * shapedic['M'] * shapedic['P'] * shapedic['Z']


def output_dict(shapedic, disagg_outputs):
    N, M, P, Z = shapedic['N'], shapedic['M'], shapedic['P'], shapedic['Z']
    dic = {}
    for out in disagg_outputs:
        shp = tuple(shapedic[key] for key in out.lower().split('_'))
        dic[out] = numpy.zeros((N, M, P) + shp + (Z,))
    return dic


@base.calculators.add('disaggregation')
class DisaggregationCalculator(base.HazardCalculator):
    """
    Classical PSHA disaggregation calculator
    """
    precalc = 'classical'
    accept_precalc = ['classical', 'disaggregation']

    def init(self):
        if self.N >= 32768:
            raise ValueError('You can disaggregate at max 32,768 sites')
        few = self.oqparam.max_sites_disagg
        if self.N > few:
            raise ValueError(
                'The number of sites is to disaggregate is %d, but you have '
                'max_sites_disagg=%d' % (self.N, few))
        super().init()

    def execute(self):
        """Performs the disaggregation"""
        return self.full_disaggregation()

    def get_curve(self, sid, rlzs):
        """
        Get the hazard curves for the given site ID and realizations.

        :param sid: site ID
        :param rlzs: a matrix of indices of shape Z
        :returns: a list of Z arrays of PoEs
        """
        poes = []
        for rlz in rlzs:
            pmap = self.pgetter.get(rlz)
            poes.append(pmap[sid].convert(self.oqparam.imtls)
                        if sid in pmap else None)
        return poes

    def check_poes_disagg(self, curves, rlzs):
        """
        Raise an error if the given poes_disagg are too small compared to
        the hazard curves.
        """
        oq = self.oqparam
        # there may be sites where the sources are too small to produce
        # an effect at the given poes_disagg
        ok_sites = []
        for sid in self.sitecol.sids:
            if all(curve is None for curve in curves[sid]):
                ok_sites.append(sid)
                continue
            bad = _check_curves(sid, rlzs[sid], curves[sid],
                                oq.imtls, oq.poes_disagg)
            if not bad:
                ok_sites.append(sid)
        if len(ok_sites) == 0:
            raise SystemExit('Cannot do any disaggregation')
        elif len(ok_sites) < self.N:
            logging.warning('Doing the disaggregation on %s', self.sitecol)
        return ok_sites

    def full_disaggregation(self):
        """
        Run the disaggregation phase.
        """
        oq = self.oqparam
        mags_by_trt = self.datastore['source_mags']
        all_edges, self.shapedic = disagg.get_edges_shapedic(
            oq, self.sitecol, mags_by_trt)
        *self.bin_edges, self.trts = all_edges
        if hasattr(self, 'csm'):
            for sg in self.csm.src_groups:
                if sg.atomic:
                    raise NotImplementedError(
                        'Atomic groups are not supported yet')
        elif self.datastore['source_info'].attrs['atomic']:
            raise NotImplementedError(
                'Atomic groups are not supported yet')

        self.full_lt = self.datastore['full_lt']
        self.poes_disagg = oq.poes_disagg or (None,)
        self.imts = list(oq.imtls)
        self.M = len(self.imts)
        ws = [rlz.weight for rlz in self.full_lt.get_realizations()]
        self.pgetter = getters.PmapGetter(
            self.datastore, ws, self.sitecol.sids)

        # build array rlzs (N, Z)
        if oq.rlz_index is None:
            Z = oq.num_rlzs_disagg or 1
            rlzs = numpy.zeros((self.N, Z), int)
            if self.R > 1:
                for sid in self.sitecol.sids:
                    curves = numpy.array(
                        [pc.array for pc in self.pgetter.get_pcurves(sid)])
                    mean = getters.build_stat_curve(
                        curves, oq.imtls, stats.mean_curve, ws)
                    rlzs[sid] = util.closest_to_ref(curves, mean.array)[:Z]
            self.datastore['best_rlzs'] = rlzs
        else:
            Z = len(oq.rlz_index)
            rlzs = numpy.zeros((self.N, Z), int)
            for z in range(Z):
                rlzs[:, z] = oq.rlz_index[z]
            self.datastore['best_rlzs'] = rlzs
        assert Z <= self.R, (Z, self.R)
        self.Z = Z
        self.rlzs = rlzs

        if oq.iml_disagg:
            # no hazard curves are needed
            self.poe_id = {None: 0}
            curves = [[None for z in range(Z)] for s in range(self.N)]
            self.ok_sites = set(self.sitecol.sids)
        else:
            self.poe_id = {poe: i for i, poe in enumerate(oq.poes_disagg)}
            curves = [self.get_curve(sid, rlzs[sid])
                      for sid in self.sitecol.sids]
            self.ok_sites = set(self.check_poes_disagg(curves, rlzs))
        self.iml4 = _iml4(rlzs, oq.iml_disagg, oq.imtls,
                          self.poes_disagg, curves)
        self.datastore['iml4'] = self.iml4
        self.datastore['poe4'] = numpy.zeros_like(self.iml4.array)

        self.save_bin_edges()
        tot = get_outputs_size(self.shapedic, oq.disagg_outputs)
        logging.info('Total output size: %s', humansize(sum(tot.values())))
        self.imldic = {}  # sid, rlz, poe, imt -> iml
        for s in self.sitecol.sids:
            iml3 = self.iml4[s]
            for z, rlz in enumerate(rlzs[s]):
                for p, poe in enumerate(self.poes_disagg):
                    for m, imt in enumerate(oq.imtls):
                        self.imldic[s, rlz, poe, imt] = iml3[m, p, z]
        return self.compute()

    def compute(self):
        """
        Submit disaggregation tasks and return the results
        """
        oq = self.oqparam
        dstore = (self.datastore.parent if self.datastore.parent
                  else self.datastore)
        mags = set()
        for trt, dset in self.datastore['source_mags'].items():
            mags.update(dset[:])
        totrups = 0
        totweight = 0
        maxsids = 0
        ns = {}  # mag -> num_sids
        for name, dset in dstore.items():
            if name.startswith('mag_'):
                nsids = numpy.array([len(sids) for sids in dset['sids_']])
                ns[name] = nsids
                maxsids = max(maxsids, nsids.max())
                totweight += nsids.sum()
                totrups += len(nsids)

        grp_ids = dstore['grp_ids'][:]
        maxw = min(int(numpy.ceil(totweight / (oq.concurrent_tasks or 1))),
                   oq.ruptures_per_block * 100)  # at maximum 50000
        rlzs_by_gsim = self.full_lt.get_rlzs_by_gsim_list(grp_ids)
        num_eff_rlzs = len(self.full_lt.sm_rlzs)
        task_inputs = []
        U, G = 0, 0
        allargs = []
        for mag in sorted(mags):
            rctx = dstore['mag_%s/rctx' % mag][:]
            for gidx, gids in enumerate(grp_ids):
                mask = rctx['gidx'] == gidx
                if mask.sum() == 0:
                    continue
                array = rctx[mask]
                nsids = ns['mag_%s' % mag][mask]
                proxies = [RuptureProxy(r, n) for r, n in zip(array,  nsids)]
                trti = gids[0] // num_eff_rlzs
                trt = self.trts[trti]
                cmaker = ContextMaker(
                    trt, rlzs_by_gsim[gidx],
                    {'truncation_level': oq.truncation_level,
                     'maximum_distance': oq.maximum_distance,
                     'collapse_level': oq.collapse_level,
                     'imtls': oq.imtls})
                G = max(G, len(cmaker.gsims))
                for blk in block_splitter(proxies, maxw,
                                          lambda p: numpy.sqrt(p.nsites)):
                    nr = len(blk)
                    U = max(U, nr)
                    arr = numpy.array([p.rec for p in blk])
                    allargs.append((dstore, arr, cmaker,
                                    self.iml4, trti, mag, self.bin_edges, oq))
                    task_inputs.append((trti, mag, nr))

        nbytes, msg = get_array_nbytes(dict(S=maxsids, M=self.M, G=G, U=U))
        logging.info('Maximum mean_std per task:\n%s', msg)

        s = self.shapedic
        size = s['dist'] * s['eps'] + s['lon'] * s['lat']
        sd = dict(N=s['N'], M=s['M'], P=s['P'], Z=s['Z'], size=size)
        sd['tasks'] = numpy.ceil(len(allargs))
        nbytes, msg = get_array_nbytes(sd)
        if nbytes > oq.max_data_transfer:
            raise ValueError(
                'Estimated data transfer too big\n%s > max_data_transfer=%s' %
                (msg, humansize(oq.max_data_transfer)))
        logging.info('Estimated data transfer:\n%s', msg)

        sd.pop('tasks')
        sd['mags_trt'] = sum(len(mags) for mags in
                             self.datastore['source_mags'].values())
        nbytes, msg = get_array_nbytes(sd)
        logging.info('Estimated memory on the master:\n%s', msg)

        dt = numpy.dtype([('trti', U8), ('mag', '|S4'), ('nrups', U32)])
        self.datastore['disagg_task'] = numpy.array(task_inputs, dt)
        self.datastore.swmr_on()
        smap = parallel.Starmap(
            compute_disagg, allargs, h5=self.datastore.hdf5)
        results = smap.reduce(self.agg_result, AccumDict(accum={}))
        return results  # imti, sid -> trti, magi -> 6D array

    def agg_result(self, acc, result):
        """
        Collect the results coming from compute_disagg into self.results.

        :param acc: dictionary sid -> trti, magi -> 6D array
        :param result: dictionary with the result coming from a task
        """
        # 7D array of shape (#distbins, #lonbins, #latbins, #epsbins, M, P, Z)
        with self.monitor('aggregating disagg matrices'):
            trti = result.pop('trti')
            magi = result.pop('magi')
            for (s, m), out in result.items():
                for k in (0, 1):
                    x = acc[s, m, k].get((trti, magi), 0)
                    acc[s, m, k][trti, magi] = agg_probs(x, out[k])
        return acc

    def save_bin_edges(self):
        """
        Save disagg-bins
        """
        b = self.bin_edges
        T = len(self.trts)
        shape = [len(bin) - 1 for bin in
                 (b[0], b[1], b[2][0], b[3][0], b[4])] + [T]
        matrix_size = numpy.prod(shape)  # 6D
        if matrix_size > 1E6:
            raise ValueError(
                'The disaggregation matrix is too large '
                '(%d elements): fix the binning!' % matrix_size)

        def a(bin_no):
            # lon/lat edges for the sites, bin_no can be 2 or 3
            num_edges = len(b[bin_no][0])
            arr = numpy.zeros((self.N, num_edges))
            for sid, edges in b[bin_no].items():
                arr[sid] = edges
            return arr
        self.datastore['disagg-bins/Mag'] = b[0]
        self.datastore['disagg-bins/Dist'] = b[1]
        self.datastore['disagg-bins/Lon'] = a(2)
        self.datastore['disagg-bins/Lat'] = a(3)
        self.datastore['disagg-bins/Eps'] = b[4]
        self.datastore['disagg-bins/TRT'] = encode(self.trts)

    def post_execute(self, results):
        """
        Save all the results of the disaggregation. NB: the number of results
        to save is #sites * #rlzs * #disagg_poes * #IMTs.

        :param results:
            a dictionary sid, imti, kind -> trti -> disagg matrix
        """
        # the DEBUG dictionary is populated only for OQ_DISTRIBUTE=no
        for sid, pnes in disagg.DEBUG.items():
            print('site %d, mean pnes=%s' % (sid, pnes))
        T = len(self.trts)
        Ma = len(self.bin_edges[0]) - 1  # num_mag_bins
        # build a dictionary s, m, k -> matrices
        results = {smk: _matrix(dic, T, Ma) for smk, dic in results.items()}
        # get the number of outputs
        shp = (self.N, len(self.poes_disagg), len(self.imts), self.Z)
        logging.info('Extracting and saving the PMFs for %d outputs '
                     '(N=%s, P=%d, M=%d, Z=%d)', numpy.prod(shp), *shp)
        with self.monitor('saving disagg results'):
            self.save_disagg_results(results)

    def save_disagg_results(self, results):
        """
        Save the computed PMFs in the datastore

        :param results:
            a dict s, m, k -> 6D-matrix of shape (T, Ma, Lo, La, P, Z) or
            (T, Ma, D, E, P, Z) depending if k is 0 or k is 1
        """
        outputs = self.oqparam.disagg_outputs
        out = output_dict(self.shapedic, outputs)
        count = numpy.zeros(len(self.sitecol), U16)
        _disagg_trt = numpy.zeros(self.N, [(trt, float) for trt in self.trts])
        for (s, m, k), mat6 in sorted(results.items()):
            imt = self.imts[m]
            for p, poe in enumerate(self.poes_disagg):
                mat5 = mat6[..., p, :]
                if k == 0 and m == 0 and poe == self.poes_disagg[-1]:
                    # mat5 has shape (T, Ma, D, E, Z)
                    _disagg_trt[s] = tuple(pprod(mat5[..., 0], axis=(1, 2, 3)))
                poe2 = pprod(mat5, axis=(0, 1, 2, 3))
                self.datastore['poe4'][s, m, p] = poe2  # shape Z
                poe_agg = poe2.mean()
                if (poe and abs(1 - poe_agg / poe) > .1 and not count[s]
                        and s in self.ok_sites):
                    logging.warning(
                        'Site #%d, IMT=%s: poe_agg=%s is quite different from '
                        'the expected poe=%s; perhaps the number of intensity '
                        'measure levels is too small?', s, imt, poe_agg, poe)
                    count[s] += 1
                mat4 = agg_probs(*mat5)  # shape (Ma D E Z) or (Ma Lo La Z)
                for key in outputs:
                    if key == 'Mag' and k == 0:
                        out[key][s, m, p, :] = pprod(mat4, axis=(1, 2))
                    elif key == 'Dist' and k == 0:
                        out[key][s, m, p, :] = pprod(mat4, axis=(0, 2))
                    elif key == 'TRT' and k == 0:
                        out[key][s, m, p, :] = pprod(mat5, axis=(1, 2, 3))
                    elif key == 'Mag_Dist' and k == 0:
                        out[key][s, m, p, :] = pprod(mat4, axis=2)
                    elif key == 'Mag_Dist_Eps' and k == 0:
                        out[key][s, m, p, :] = mat4
                    elif key == 'Lon_Lat' and k == 1:
                        out[key][s, m, p, :] = pprod(mat4, axis=0)
                    elif key == 'Mag_Lon_Lat' and k == 1:
                        out[key][s, m, p, :] = mat4
                    elif key == 'Lon_Lat_TRT' and k == 1:
                        out[key][s, m, p, :] = pprod(mat5, axis=1).transpose(
                            1, 2, 0, 3)  # T Lo La Z -> Lo La T Z
                    # shape NMP..Z
        self.datastore['disagg'] = out
        # below a dataset useful for debugging, at minimum IMT and maximum RP
        self.datastore['_disagg_trt'] = _disagg_trt
